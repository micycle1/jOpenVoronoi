package org.rogach.jopenvoronoi.pocket;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;
import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Fills a closed 2D shape with a single continuous spiral polyline.
 *
 * <p>
 * The spiral starts at the widest interior point of the shape (the center of
 * its largest inscribed circle), winds outward in evenly spaced laps, and ends
 * exactly on the shape's outline. Neighboring laps are never further apart
 * than the configured {@linkplain #setStepOver(double) step-over}, so the
 * spiral visually "fills" the shape — dense hatching for a small step-over,
 * airy rings for a large one. Concave shapes, lobes and spikes are all
 * followed; each connected region of the shape's medial axis yields one
 * spiral.
 *
 * <h2>Typical usage</h2>
 *
 * <pre>{@code
 * // 1. build a Voronoi diagram of the closed outline (CCW polygon)
 * VoronoiDiagram vd = new VoronoiDiagram();
 * List<Vertex> vs = new ArrayList<>();
 * for (Point p : outline) {
 * 	vs.add(vd.insertPointSite(p));
 * }
 * for (int i = 0; i < vs.size(); i++) {
 * 	vd.insertLineSite(vs.get(i), vs.get((i + 1) % vs.size()));
 * }
 *
 * // 2. reduce it to the interior medial axis (required!)
 * vd.filter(new PolygonInteriorFilter(true));
 * vd.filter(new MedialAxisFilter());
 *
 * // 3. spiral
 * SpiralToolPath spiral = new SpiralToolPath(vd.getDiagram());
 * spiral.setStepOver(0.01); // lap spacing, in input coordinates
 * spiral.run();
 * List<List<Point>> paths = spiral.getToolPathComponents();
 * }</pre>
 *
 * <p>
 * The constructor rejects diagrams that have not been run through
 * {@code PolygonInteriorFilter} and {@code MedialAxisFilter}. As everywhere in
 * jOpenVoronoi, input coordinates should live roughly within the unit circle.
 *
 * <h2>Knobs</h2>
 * <ul>
 * <li>{@link #setStepOver(double)} — lap spacing; the dominant visual
 * control.</li>
 * <li>{@link #setSmoothing(double)} — how much of the step-over budget is
 * spent rounding off corners (default {@code 0.5}); {@code 0} yields the raw
 * angular spiral.</li>
 * <li>{@link #setSamplingLength(double)} — geometric fidelity of the
 * underlying medial-axis discretization; rarely needs changing.</li>
 * <li>{@link #setStoreLaps(boolean)} / {@link #setStoreWavefronts(boolean)} —
 * retain intermediate geometry (concentric rings), interesting both for
 * debugging and as artwork in its own right.</li>
 * </ul>
 *
 * <h2>Algorithm</h2>
 *
 * Implements the linearized medial-axis spiral of Held &amp; de Lorenzo
 * (EuroCG 2017) / Held &amp; Spielberger (CAD 2009), following the Spielberger
 * thesis: the medial axis is discretized adaptively (both along the axis and
 * by the spacing of clearance-line endpoints on the boundary, which keeps
 * concave regions faithful), an impulse propagated from the root defines
 * wavefronts, laps are interpolated between wavefronts and chained into one
 * spiral, and the raw path is finally faired inside a tolerance band whose
 * width is taken out of the step-over budget — so the smoothed path still
 * respects the configured step-over and never leaves the shape.
 *
 * <p>
 * Intended for simple pockets (tree-like medial axis). Islands / cyclic
 * medial-axis components are not supported.
 */
public class SpiralToolPath {

	private static final double EPS = 1e-9;

	/** Maximum recursive refinement depth per seed interval of an MA edge. */
	private static final int MAX_REFINE_DEPTH = 8;
	/** Maximum number of uniform seed intervals per MA edge. */
	private static final int MAX_SEED_SEGMENTS = 512;
	/** Fraction of the local clearance that smoothing is allowed to consume. */
	private static final double CLEARANCE_SAFETY = 0.9;
	/** Upper bound on densified point count fed to the smoothing relaxation. */
	private static final int MAX_SMOOTHING_POINTS = 500_000;
	private static final int SMOOTHING_ITERATIONS = 200;
	private static final double SMOOTHING_RELAX = 0.5;

	private final HalfEdgeDiagram g;
	private final List<Edge> maEdges = new ArrayList<>();

	private double maxStepOver = 0.05;
	private double smoothing = 0.5;
	private double samplingLength = 0.0; // 0 = auto (effective step-over)
	private int samplesPerCurveEdge = 10;
	private int maxLaps = 2000;

	private boolean storeWavefronts = false;
	private boolean storeLaps = false;

	private final List<List<Point>> toolPathComponents = new ArrayList<>();
	private final List<List<Point>> rawToolPathComponents = new ArrayList<>();
	private final List<List<List<Point>>> lapComponents = new ArrayList<>();
	private final List<List<List<Point>>> wavefrontComponents = new ArrayList<>();

	/**
	 * Creates a spiral generator for the given filtered Voronoi diagram.
	 * <p>
	 * The diagram must already be reduced to the interior medial axis of the
	 * shape, i.e. both {@code vd.filter(new PolygonInteriorFilter(true))} and
	 * {@code vd.filter(new MedialAxisFilter())} must have been applied — see the
	 * class example. Anything else would spiral over the whole plane rather than
	 * the inside of the shape, so it is rejected eagerly.
	 *
	 * @param g the half-edge diagram obtained from
	 *          {@code VoronoiDiagram.getDiagram()} after filtering
	 * @throws IllegalStateException if the diagram still contains edges that the
	 *                               two required filters would have removed
	 */
	public SpiralToolPath(HalfEdgeDiagram g) {
		this.g = g;

		boolean unfilteredInterior = false;
		boolean unfilteredExterior = false;
		for (Edge e : g.edges) {
			if (!e.valid) {
				continue;
			}
			if (e.type == EdgeType.OUTEDGE) {
				unfilteredExterior = true;
			} else if (e.type == EdgeType.SEPARATOR) {
				unfilteredInterior = true;
			} else if (isMedialAxisEdge(e)) {
				maEdges.add(e);
			}
		}
		if (unfilteredExterior || unfilteredInterior) {
			throw new IllegalStateException("Diagram has not been reduced to its interior medial axis"
					+ (unfilteredExterior ? " (exterior edges are still valid; PolygonInteriorFilter missing?)" : "")
					+ (unfilteredInterior ? " (separator edges are still valid; MedialAxisFilter missing?)" : "")
					+ ". Apply vd.filter(new PolygonInteriorFilter(true)) followed by vd.filter(new MedialAxisFilter()) before building a spiral.");
		}
	}

	/**
	 * Maximum distance between neighboring laps of the spiral, in the same units
	 * as the input coordinates. This is the main visual control: smaller values
	 * give a denser spiral (and proportionally more laps and points). Default
	 * {@code 0.05}.
	 *
	 * @param stepOver lap spacing, {@code > 0}
	 */
	public void setStepOver(double stepOver) {
		if (stepOver <= 0) {
			throw new IllegalArgumentException("stepOver must be > 0");
		}
		this.maxStepOver = stepOver;
	}

	/**
	 * Fraction of the step-over budget spent on rounding off corners, in
	 * {@code [0, 0.8]}. Default {@code 0.5}.
	 * <p>
	 * The raw spiral is a polyline with visible kinks wherever the medial axis
	 * branches. Smoothing generates the spiral with a tighter effective step-over
	 * of {@code (1 - smoothing) * stepOver} and then fairs it, letting every point
	 * drift up to {@code smoothing * stepOver / 2} from its raw position (never
	 * across the shape's outline). The combined result still respects
	 * {@link #setStepOver(double)} everywhere.
	 * <p>
	 * Higher values look calmer but produce more laps ({@code 0.5} doubles the lap
	 * count relative to {@code 0}). {@code 0} disables smoothing entirely, which
	 * can be attractive when the faceted look is wanted — the raw spiral is also
	 * always available from {@link #getRawToolPathComponents()}.
	 *
	 * @param smoothing smoothing fraction in {@code [0, 0.8]}
	 */
	public void setSmoothing(double smoothing) {
		if (smoothing < 0 || smoothing > 0.8) {
			throw new IllegalArgumentException("smoothing must be in [0, 0.8]");
		}
		this.smoothing = smoothing;
	}

	/**
	 * Sampling length used to discretize the medial axis: sample nodes are placed
	 * so that neither the distance between consecutive nodes nor the distance
	 * between their contact points on the shape's outline exceeds this value. The
	 * second criterion is what keeps the spiral faithful inside concave features.
	 * <p>
	 * By default ({@code 0}) the sampling length tracks the effective step-over,
	 * which is almost always right: fidelity automatically increases with spiral
	 * density. Set it explicitly only to trade accuracy for speed (larger values)
	 * or to sharpen very fine boundary detail (smaller values, at a roughly linear
	 * cost in time and memory).
	 *
	 * @param samplingLength sampling length in input units, or {@code 0} for
	 *                       automatic
	 */
	public void setSamplingLength(double samplingLength) {
		if (samplingLength < 0) {
			throw new IllegalArgumentException("samplingLength must be >= 0");
		}
		this.samplingLength = samplingLength;
	}

	/**
	 * Minimum number of uniform seed samples per curved medial-axis edge, refined
	 * adaptively afterwards. Straight edges are seeded from their length and the
	 * sampling length directly. Default {@code 10}; rarely needs changing since
	 * the adaptive refinement dominates the result.
	 *
	 * @param n seed sample count, {@code >= 1}
	 */
	public void setSamplesPerCurveEdge(int n) {
		if (n < 1) {
			throw new IllegalArgumentException("samplesPerCurveEdge must be >= 1");
		}
		this.samplesPerCurveEdge = n;
	}

	/**
	 * Safety valve: {@link #run()} refuses to build more laps than this per
	 * component (default {@code 2000}), which guards against an accidentally tiny
	 * step-over on a large shape. Raise it deliberately for genuinely dense work.
	 *
	 * @param maxLaps maximum lap count, {@code >= 1}
	 */
	public void setMaxLaps(int maxLaps) {
		if (maxLaps < 1) {
			throw new IllegalArgumentException("maxLaps must be >= 1");
		}
		this.maxLaps = maxLaps;
	}

	/**
	 * When enabled, {@link #run()} also retains the sampled wavefronts — the
	 * closed rings the spiral interpolates between — retrievable from
	 * {@link #getWavefrontComponents()}. Off by default (costs memory).
	 */
	public void setStoreWavefronts(boolean storeWavefronts) {
		this.storeWavefronts = storeWavefronts;
	}

	/**
	 * When enabled, {@link #run()} also retains the individual laps of the raw
	 * spiral, retrievable from {@link #getLapComponents()}. Off by default (costs
	 * memory).
	 */
	public void setStoreLaps(boolean storeLaps) {
		this.storeLaps = storeLaps;
	}

	/**
	 * Returns the finished spiral, one polyline per connected medial-axis
	 * component of the shape.
	 * <p>
	 * Each polyline starts at the widest interior point, winds outward and ends on
	 * the shape's outline; consecutive points are close together, so it can be
	 * drawn directly as a line strip. If smoothing is enabled (default) this is
	 * the faired path; the angular raw path is available from
	 * {@link #getRawToolPathComponents()}.
	 * <p>
	 * Populated by {@link #run()}.
	 *
	 * @return one spiral polyline per component; empty before {@link #run()}
	 */
	public List<List<Point>> getToolPathComponents() {
		return toolPathComponents;
	}

	/**
	 * Returns the raw (unsmoothed) spiral, one polyline per connected medial-axis
	 * component — same geometry the smoothing stage starts from, with its faceted,
	 * angular character intact. Populated by {@link #run()} regardless of the
	 * smoothing setting.
	 *
	 * @return one raw spiral polyline per component; empty before {@link #run()}
	 */
	public List<List<Point>> getRawToolPathComponents() {
		return rawToolPathComponents;
	}

	/**
	 * Returns the individual laps of the raw spiral for each component. A lap is
	 * one full outward turn; lap {@code 0} is the innermost and the last lap runs
	 * along the shape's outline. The nesting is
	 * {@code component -> lap -> polyline points}.
	 * <p>
	 * Only retained when {@link #setStoreLaps(boolean)} was enabled before
	 * {@link #run()}; otherwise this returns an empty list.
	 *
	 * @return lap polylines per component, or an empty list when lap storage is
	 *         disabled
	 */
	public List<List<List<Point>>> getLapComponents() {
		return lapComponents;
	}

	/**
	 * Returns the sampled wavefronts for each component: closed concentric rings
	 * marching from the interior to the outline, evaluated at evenly spaced times
	 * of the impulse-propagation model. The nesting is
	 * {@code component -> wavefront -> polyline points}, ordered inside-out.
	 * <p>
	 * Only retained when {@link #setStoreWavefronts(boolean)} was enabled before
	 * {@link #run()}; otherwise this returns an empty list.
	 *
	 * @return wavefront polylines per component, or an empty list when wavefront
	 *         storage is disabled
	 */
	public List<List<List<Point>>> getWavefrontComponents() {
		return wavefrontComponents;
	}

	/**
	 * Computes the spiral(s). Results are exposed through
	 * {@link #getToolPathComponents()} (and the other getters); calling this again
	 * after changing settings recomputes everything from scratch.
	 *
	 * @throws IllegalArgumentException if the configured step-over would require
	 *                                  more than {@link #setMaxLaps(int)} laps
	 */
	public void run() {
		toolPathComponents.clear();
		rawToolPathComponents.clear();
		lapComponents.clear();
		wavefrontComponents.clear();

		double effectiveStepOver = maxStepOver * (1.0 - smoothing);
		double lambda = Math.max((samplingLength > 0) ? samplingLength : effectiveStepOver, 1e-6);

		List<List<Edge>> components = findMedialAxisComponents();
		for (List<Edge> component : components) {
			if (component.isEmpty()) {
				continue;
			}

			Set<Edge> componentCanonicalSet = identitySet();
			componentCanonicalSet.addAll(component);

			DiscreteTree discreteTree = buildDiscreteTree(component, componentCanonicalSet, lambda);
			if (discreteTree == null || discreteTree.root == null) {
				continue;
			}

			sortChildrenByEmbedding(discreteTree.root);
			computeHeightsAndStartTimes(discreteTree.root);
			buildAncestorTable(discreteTree.root);

			List<DiscreteNode> leafOrder = new ArrayList<>();
			collectLeaves(discreteTree.root, leafOrder);
			if (leafOrder.isEmpty()) {
				continue;
			}

			List<Branch> branches = buildBranches(leafOrder, discreteTree.root);

			int m = Math.max(1, (int) Math.ceil(discreteTree.root.height / effectiveStepOver));
			if (m > maxLaps) {
				throw new IllegalArgumentException("Requested step-over creates too many laps: " + m + " (height=" + discreteTree.root.height
						+ ", effectiveStepOver=" + effectiveStepOver + ")");
			}
			double dt = 1.0 / m;

			if (storeWavefronts) {
				List<List<TreePosition>> wavefronts = new ArrayList<>();
				for (int i = 0; i <= m; i++) {
					wavefronts.add(evaluateWavefront(branches, i * dt));
				}
				wavefrontComponents.add(toPointChains(wavefronts));
			}

			List<List<TreePosition>> laps = buildLaps(branches, m, dt);

			if (storeLaps) {
				lapComponents.add(toPointChains(laps));
			}

			List<TreePosition> connected = connectLaps(laps);
			rawToolPathComponents.add(toPoints(connected));
			if (smoothing > 0) {
				toolPathComponents.add(smoothPath(connected));
			} else {
				toolPathComponents.add(toPoints(connected));
			}
		}
	}

	// -------------------------------------------------------------------------
	// Build rooted original MA tree
	// -------------------------------------------------------------------------

	private Vertex selectRootVertex(List<Edge> component) {
		Vertex best = null;
		double bestClr = -1;

		Set<Vertex> seen = identitySet();
		for (Edge e : component) {
			if (seen.add(e.source) && e.source.dist() > bestClr) {
				bestClr = e.source.dist();
				best = e.source;
			}
			if (seen.add(e.target) && e.target.dist() > bestClr) {
				bestClr = e.target.dist();
				best = e.target;
			}
		}
		return best;
	}

	// -------------------------------------------------------------------------
	// Build discrete tree
	// -------------------------------------------------------------------------

	private DiscreteTree buildDiscreteTree(List<Edge> component, Set<Edge> componentCanonicalSet, double lambda) {
		Vertex rootVertex = selectRootVertex(component);
		if (rootVertex == null) {
			return null;
		}

		DiscreteTree tree = new DiscreteTree();
		tree.root = createVertexNode(rootVertex, componentCanonicalSet, tree);

		Set<Vertex> visited = identitySet();
		visited.add(rootVertex);

		Deque<VertexFrame> stack = new ArrayDeque<>();
		stack.push(new VertexFrame(rootVertex, tree.root));

		while (!stack.isEmpty()) {
			VertexFrame frame = stack.pop();
			Vertex v = frame.vertex;
			DiscreteNode discreteParent = frame.node;

			for (Edge oe : v.outEdges) {
				Edge canonical = canonicalRepresentative(oe, componentCanonicalSet);
				if (canonical == null) {
					continue;
				}

				Vertex neighbor = otherEndpoint(canonical, v);
				if (neighbor == null || visited.contains(neighbor)) {
					continue;
				}

				Edge oriented = orientFrom(canonical, v);
				if (oriented == null) {
					continue;
				}

				DiscreteNode prev = discreteParent;
				for (double u : interiorSampleParams(oriented, lambda)) {
					DiscreteNode sampleNode = createSampleNode(oriented, u, tree);
					prev = connectOrMerge(prev, sampleNode, tree);
				}

				DiscreteNode childVertexNode = createVertexNode(neighbor, componentCanonicalSet, tree);
				// connectOrMerge may merge childVertexNode into prev (coincident
				// vertices, e.g. multi-way junctions built from zero-length edges);
				// the surviving node must carry the traversal or the subtree behind
				// the junction ends up disconnected from the root.
				prev = connectOrMerge(prev, childVertexNode, tree);

				visited.add(neighbor);
				stack.push(new VertexFrame(neighbor, prev));
			}
		}

		attachClearanceLeaves(tree);
		return tree;
	}

	// -------------------------------------------------------------------------
	// Adaptive MA-edge discretization (thesis section 3.5)
	// -------------------------------------------------------------------------

	private boolean isSampleable(EdgeType type) {
		return type == EdgeType.LINE || type == EdgeType.LINELINE || type == EdgeType.PARA_LINELINE || type == EdgeType.PARABOLA;
	}

	/**
	 * Returns interior parameters (strictly between 0 and 1, ascending) at which
	 * the given oriented MA edge should be sampled. The edge is first seeded
	 * uniformly from its chord length (curved edges additionally honour
	 * {@link #setSamplesPerCurveEdge(int)}), then every interval is refined
	 * recursively while either the MA chord or the spacing of the boundary touch
	 * points of its endpoints exceeds {@code lambda}. The touch-point criterion is
	 * what keeps the discretization dense around concave boundary features.
	 */
	private List<Double> interiorSampleParams(Edge e, double lambda) {
		if (!isSampleable(e.type)) {
			return Collections.emptyList();
		}

		double chord = distance(e.source.position, e.target.position);
		int seed = (int) Math.ceil(chord / lambda);
		if (e.type == EdgeType.PARABOLA) {
			seed = Math.max(seed, samplesPerCurveEdge);
		}
		seed = Math.max(1, Math.min(seed, MAX_SEED_SEGMENTS));

		List<Double> out = new ArrayList<>();
		EdgeSample prev = edgeSample(e, 0.0);
		for (int i = 1; i <= seed; i++) {
			double u = i / (double) seed;
			EdgeSample cur = edgeSample(e, u);
			refine(e, prev, cur, lambda, 0, out);
			if (i < seed) {
				out.add(u);
			}
			prev = cur;
		}
		return out;
	}

	private void refine(Edge e, EdgeSample s0, EdgeSample s1, double lambda, int depth, List<Double> out) {
		if (depth >= MAX_REFINE_DEPTH || !needsSplit(s0, s1, lambda)) {
			return;
		}
		double um = 0.5 * (s0.u + s1.u);
		EdgeSample mid = edgeSample(e, um);
		refine(e, s0, mid, lambda, depth + 1, out);
		out.add(um);
		refine(e, mid, s1, lambda, depth + 1, out);
	}

	private boolean needsSplit(EdgeSample s0, EdgeSample s1, double lambda) {
		if (distance(s0.point, s1.point) > lambda) {
			return true;
		}
		if (distance(s0.touchA, s1.touchA) > lambda) {
			return true;
		}
		return distance(s0.touchB, s1.touchB) > lambda;
	}

	private EdgeSample edgeSample(Edge e, double u) {
		EdgeSample s = new EdgeSample();
		s.u = u;
		s.point = e.micSample(u).getKey();
		Point[] tp = e.micTouchPoints(u);
		s.touchA = tp[0];
		s.touchB = tp[1];
		return s;
	}

	private static final class EdgeSample {
		double u;
		Point point;
		Point touchA;
		Point touchB;
	}

	private static final class VertexFrame {
		final Vertex vertex;
		final DiscreteNode node;

		VertexFrame(Vertex vertex, DiscreteNode node) {
			this.vertex = vertex;
			this.node = node;
		}
	}

	private DiscreteNode createVertexNode(Vertex v, Set<Edge> componentCanonicalSet, DiscreteTree tree) {
		DiscreteNode n = new DiscreteNode();
		n.point = v.position;
		n.clearance = v.dist();

		for (Edge oe : v.outEdges) {
			Edge canonical = canonicalRepresentative(oe, componentCanonicalSet);
			if (canonical == null || !isSampleable(oe.type)) {
				continue;
			}
			addTouchPoints(n.touchPoints, oe.micTouchPoints(0.0));
		}

		tree.maNodes.add(n);
		return n;
	}

	private DiscreteNode createSampleNode(Edge e, double u, DiscreteTree tree) {
		Entry<Point, Double> sample = e.micSample(u);

		DiscreteNode n = new DiscreteNode();
		n.point = sample.getKey();
		n.clearance = sample.getValue();
		addTouchPoints(n.touchPoints, e.micTouchPoints(u));

		tree.maNodes.add(n);
		return n;
	}

	private DiscreteNode connectOrMerge(DiscreteNode parent, DiscreteNode child, DiscreteTree tree) {
		double len = distance(parent.point, child.point);

		if (len <= EPS) {
			addTouchPoints(parent.touchPoints, child.touchPoints.toArray(new Point[0]));
			// child is always the most recently created node
			int last = tree.maNodes.size() - 1;
			if (last >= 0 && tree.maNodes.get(last) == child) {
				tree.maNodes.remove(last);
			} else {
				tree.maNodes.remove(child);
			}
			return parent;
		}

		DiscreteEdge de = new DiscreteEdge();
		de.parent = parent;
		de.child = child;
		de.length = len;

		child.parent = parent;
		child.parentEdge = de;
		parent.children.add(de);

		return child;
	}

	private void attachClearanceLeaves(DiscreteTree tree) {
		for (DiscreteNode center : tree.maNodes) {
			for (Point tp : center.touchPoints) {
				if (tp == null || distance(center.point, tp) <= EPS) {
					continue;
				}

				DiscreteNode leaf = new DiscreteNode();
				leaf.point = tp;
				leaf.clearance = 0.0;
				leaf.boundaryLeaf = true;

				DiscreteEdge edge = new DiscreteEdge();
				edge.parent = center;
				edge.child = leaf;
				edge.length = distance(center.point, leaf.point);

				leaf.parent = center;
				leaf.parentEdge = edge;
				center.children.add(edge);
			}
		}
	}

	// -------------------------------------------------------------------------
	// Tree ordering / height / times
	// -------------------------------------------------------------------------

	private void sortChildrenByEmbedding(DiscreteNode root) {
		Deque<DiscreteNode> stack = new ArrayDeque<>();
		stack.push(root);

		while (!stack.isEmpty()) {
			DiscreteNode node = stack.pop();

			if (node.children.size() > 1) {
				if (node.parent == null) {
					for (DiscreteEdge e : node.children) {
						e.sortKey = angleOf(node.point, e.child.point);
					}
				} else {
					final double base = angleOf(node.point, node.parent.point);
					for (DiscreteEdge e : node.children) {
						e.sortKey = ccwDelta(base, angleOf(node.point, e.child.point));
					}
				}
				node.children.sort(Comparator.comparingDouble(e -> e.sortKey));
			}

			for (DiscreteEdge e : node.children) {
				stack.push(e.child);
			}
		}
	}

	private void computeHeightsAndStartTimes(DiscreteNode root) {
		List<DiscreteNode> order = new ArrayList<>();
		Deque<DiscreteNode> stack = new ArrayDeque<>();
		stack.push(root);

		while (!stack.isEmpty()) {
			DiscreteNode n = stack.pop();
			order.add(n);
			for (DiscreteEdge e : n.children) {
				stack.push(e.child);
			}
		}

		// postorder: heights
		for (int i = order.size() - 1; i >= 0; i--) {
			DiscreteNode node = order.get(i);

			if (node.children.isEmpty()) {
				node.height = node.boundaryLeaf ? 0.0 : node.clearance;
				node.bestLeaf = node;
				node.bestChild = null;
				continue;
			}

			double best = -Double.MAX_VALUE;
			DiscreteEdge bestEdge = null;
			DiscreteNode bestLeaf = null;

			for (DiscreteEdge e : node.children) {
				double h = e.length + e.child.height;
				if (h > best) {
					best = h;
					bestEdge = e;
					bestLeaf = e.child.bestLeaf;
				}
			}

			node.height = best;
			node.bestChild = bestEdge;
			node.bestLeaf = bestLeaf;
		}

		// preorder: start times
		root.startTime = 0.0;
		for (DiscreteNode parent : order) {
			for (DiscreteEdge e : parent.children) {
				double L = e.length;
				double hq = e.child.height;
				e.child.startTime = (parent.startTime * hq + L) / (hq + L);
			}
		}
	}

	/**
	 * Assigns every node its depth and binary-lifting ancestor table
	 * ({@code up[k]} is the {@code 2^k}-th ancestor). This lets
	 * {@link #pointOnBranchAtTime} locate the active edge of any root-leaf path in
	 * O(log depth) without materializing per-branch edge arrays.
	 */
	private void buildAncestorTable(DiscreteNode root) {
		Deque<DiscreteNode> stack = new ArrayDeque<>();
		root.depth = 0;
		root.up = null;
		stack.push(root);

		while (!stack.isEmpty()) {
			DiscreteNode node = stack.pop();
			for (DiscreteEdge e : node.children) {
				DiscreteNode c = e.child;
				int d = node.depth + 1;
				c.depth = d;
				DiscreteNode[] up = new DiscreteNode[32 - Integer.numberOfLeadingZeros(d)]; // floor(log2(d)) + 1
				up[0] = node;
				for (int k = 1; k < up.length; k++) {
					up[k] = up[k - 1].up[k - 1];
				}
				c.up = up;
				stack.push(c);
			}
		}
	}

	private void collectLeaves(DiscreteNode root, List<DiscreteNode> out) {
		Deque<DiscreteNode> stack = new ArrayDeque<>();
		stack.push(root);

		while (!stack.isEmpty()) {
			DiscreteNode node = stack.pop();
			if (node.children.isEmpty()) {
				out.add(node);
				continue;
			}
			for (int i = node.children.size() - 1; i >= 0; i--) {
				stack.push(node.children.get(i).child);
			}
		}
	}

	// -------------------------------------------------------------------------
	// Branches / wavefronts / laps
	// -------------------------------------------------------------------------

	private List<List<TreePosition>> buildLaps(List<Branch> branches, int m, double dt) {
		List<List<TreePosition>> out = new ArrayList<>(m + 1);

		// First lap R0 from w(t0), w(t1)
		List<TreePosition> firstLap = buildEndpointLap(branches, 0.0, dt);

		// Last lap Rm from w(t_{m-1}), w(tm), then clipped as in section 3.11
		List<TreePosition> rawLastLap = buildEndpointLap(branches, 1.0 - dt, 1.0);
		List<TreePosition> lastLap = adjustLastLap(firstLap, rawLastLap, m, dt);

		out.add(firstLap);

		// Direct branch-indexed lookup instead of HashMap<Branch, TreePosition>
		TreePosition[] lastByBranch = new TreePosition[branches.size()];
		for (TreePosition p : lastLap) {
			lastByBranch[p.branch.index] = p;
		}

		List<TreePosition> currentLap = firstLap;

		// Build R1 ... R_{m-1}
		for (int i = 0; i <= m - 2; i++) {
			List<TreePosition> nextLap = new ArrayList<>(currentLap.size());
			int remaining = m - i;

			for (TreePosition pi : currentLap) {
				// Longest descendant branch starting at pi
				Branch longest = longestDescendantBranch(pi);

				// Corresponding point on the last lap
				TreePosition pm = (longest != null) ? lastByBranch[longest.index] : null;
				if (pm == null) {
					pm = lastByBranch[pi.branch.index];
				}
				if (pm == null) {
					pm = pi;
				}

				// ν_i(pi) = (h(pi) - h(pm)) / (m - i)
				double nu = Math.max(0.0, (pi.height - pm.height) / remaining);

				// Move outwards along the longest branch by ν_i(pi)
				TreePosition q = advanceAlongLongestBranch(pi, nu);

				// Then evaluate the same branch as pi at time t_q
				nextLap.add(pointOnBranchAtTime(pi.branch, q.time));
			}

			out.add(nextLap);
			currentLap = nextLap;
		}

		out.add(lastLap);
		return out;
	}

	private List<Branch> buildBranches(List<DiscreteNode> leafOrder, DiscreteNode root) {
		List<Branch> out = new ArrayList<>(leafOrder.size());

		for (int i = 0; i < leafOrder.size(); i++) {
			DiscreteNode leaf = leafOrder.get(i);

			Branch b = new Branch();
			b.index = i;
			b.root = root;
			b.leaf = leaf;

			leaf.branch = b;
			out.add(b);
		}

		return out;
	}

	private List<TreePosition> buildEndpointLap(List<Branch> branches, double t0, double t1) {
		List<TreePosition> wf = evaluateWavefront(branches, t1);
		int n = wf.size();

		double[] cum = new double[n];
		for (int i = 1; i < n; i++) {
			cum[i] = cum[i - 1] + distance(wf.get(i - 1).point, wf.get(i).point);
		}

		double perimeter = 0.0;
		if (n > 1) {
			perimeter = cum[n - 1] + distance(wf.get(n - 1).point, wf.get(0).point);
		}

		double span = t1 - t0;
		List<TreePosition> lap = new ArrayList<>(n);

		if (perimeter <= EPS) {
			for (Branch branch : branches) {
				lap.add(pointOnBranchAtTime(branch, t0));
			}
			return lap;
		}

		for (int i = 0; i < n; i++) {
			double tau = t0 + span * (cum[i] / perimeter);
			lap.add(pointOnBranchAtTime(branches.get(i), tau));
		}
		return lap;
	}

	private List<TreePosition> adjustLastLap(List<TreePosition> firstLap, List<TreePosition> lastLap, int m, double dt) {
		List<TreePosition> out = new ArrayList<>(lastLap.size());

		for (int i = 0; i < lastLap.size(); i++) {
			TreePosition p = firstLap.get(i);
			TreePosition q = lastLap.get(i);
			double maxTime = p.time + (m - 1) * dt;
			out.add(q.time > maxTime + EPS ? pointOnBranchAtTime(q.branch, maxTime) : q);
		}

		return out;
	}

	private List<TreePosition> evaluateWavefront(List<Branch> branches, double time) {
		List<TreePosition> wf = new ArrayList<>(branches.size());
		for (Branch b : branches) {
			wf.add(pointOnBranchAtTime(b, time));
		}
		return wf;
	}

	private Branch longestDescendantBranch(TreePosition p) {
		DiscreteNode bestLeaf = (p.node != null) ? p.node.bestLeaf : (p.edge != null) ? p.edge.child.bestLeaf : null;

		Branch b = (bestLeaf != null) ? bestLeaf.branch : null;
		return (b != null) ? b : p.branch;
	}

	private TreePosition advanceAlongLongestBranch(TreePosition start, double dist) {
		dist = Math.max(0.0, dist);
		Branch longest = longestDescendantBranch(start);

		if (start.node != null) {
			if (dist <= EPS) {
				return positionAtNode(longest, start.node);
			}

			DiscreteNode current = start.node;
			while (true) {
				DiscreteEdge next = current.bestChild;
				if (next == null) {
					return positionAtNode(longest, current);
				}
				if (dist <= next.length + EPS) {
					return positionOnEdge(longest, next, Math.min(next.length, dist),
							current.startTime + (next.child.startTime - current.startTime) * (Math.min(next.length, dist) / next.length));
				}
				dist -= next.length;
				current = next.child;
			}
		}

		DiscreteEdge e = start.edge;
		double remain = e.length - start.offset;

		if (dist <= remain + EPS) {
			double offset = Math.min(e.length, start.offset + dist);
			double alpha = (e.length <= EPS) ? 1.0 : offset / e.length;
			double time = e.parent.startTime + alpha * (e.child.startTime - e.parent.startTime);
			return positionOnEdge(longest, e, offset, time);
		}

		dist -= remain;
		DiscreteNode current = e.child;

		while (true) {
			DiscreteEdge next = current.bestChild;
			if (next == null) {
				return positionAtNode(longest, current);
			}
			if (dist <= next.length + EPS) {
				double offset = Math.min(next.length, dist);
				double alpha = (next.length <= EPS) ? 1.0 : offset / next.length;
				double time = next.parent.startTime + alpha * (next.child.startTime - next.parent.startTime);
				return positionOnEdge(longest, next, offset, time);
			}
			dist -= next.length;
			current = next.child;
		}
	}

	private TreePosition pointOnBranchAtTime(Branch branch, double time) {
		time = clamp(time, 0.0, 1.0);

		DiscreteNode leaf = branch.leaf;
		if (leaf.up == null) { // single-node component: leaf is the root
			return positionAtNode(branch, leaf);
		}

		if (time <= branch.root.startTime + EPS) {
			return positionAtNode(branch, branch.root);
		}
		if (time >= leaf.startTime - EPS) {
			return positionAtNode(branch, leaf);
		}

		// Start times increase strictly from root to leaf, so the nodes whose
		// start time exceeds `time` form a suffix of the path. Binary-lift to the
		// topmost such node; its parent edge is the active edge at `time`.
		DiscreteNode cur = leaf;
		for (int k = cur.up.length - 1; k >= 0; k--) {
			if (k < cur.up.length) {
				DiscreteNode anc = cur.up[k];
				if (anc.startTime > time + EPS) {
					cur = anc;
				}
			}
		}

		DiscreteEdge e = cur.parentEdge;
		double t0 = e.parent.startTime;
		double t1 = e.child.startTime;
		double alpha = (t1 - t0 <= EPS) ? 1.0 : clamp((time - t0) / (t1 - t0), 0.0, 1.0);

		return positionOnEdge(branch, e, alpha * e.length, time);
	}

	private TreePosition positionAtNode(Branch branch, DiscreteNode node) {
		TreePosition p = new TreePosition();
		p.branch = branch;
		p.node = node;
		p.point = node.point;
		p.time = node.startTime;
		p.height = node.height;
		return p;
	}

	private TreePosition positionOnEdge(Branch branch, DiscreteEdge edge, double offset, double time) {
		offset = clamp(offset, 0.0, edge.length);

		TreePosition p = new TreePosition();
		p.branch = branch;
		p.edge = edge;
		p.offset = offset;
		p.point = interpolate(edge.parent.point, edge.child.point, edge.length <= EPS ? 0.0 : offset / edge.length);
		p.time = time;
		p.height = edge.child.height + (edge.length - offset);
		return p;
	}

	// -------------------------------------------------------------------------
	// Smoothing (thesis chapter 9: tolerance-band fairing of the raw path)
	// -------------------------------------------------------------------------

	/**
	 * Conservative lower bound on the distance from a tree position to the pocket
	 * boundary, from the 1-Lipschitz property of the clearance function.
	 */
	private double clearanceBound(TreePosition p) {
		if (p.node != null) {
			return Math.max(0.0, p.node.clearance);
		}
		DiscreteEdge e = p.edge;
		double fromParent = e.parent.clearance - p.offset;
		double fromChild = e.child.clearance - (e.length - p.offset);
		return Math.max(0.0, Math.max(fromParent, fromChild));
	}

	/**
	 * Fairs the connected raw spiral with an iterative, constrained Laplacian
	 * relaxation. Every point may deviate from its raw position by at most
	 * {@code smoothing * stepOver / 2}, further capped by the local clearance so
	 * the smoothed path stays strictly inside the pocket; points on the boundary
	 * lap (clearance zero) therefore stay pinned.
	 */
	private List<Point> smoothPath(List<TreePosition> path) {
		int n = path.size();
		double tol = 0.5 * smoothing * maxStepOver;
		if (n < 3 || tol <= EPS) {
			return toPoints(path);
		}

		// Densification spacing: enough points to round corners within `tol`,
		// bounded so pathological inputs cannot explode.
		double totalLen = 0;
		for (int i = 1; i < n; i++) {
			totalLen += distance(path.get(i - 1).point, path.get(i).point);
		}
		double h = Math.max(0.5 * tol, totalLen / MAX_SMOOTHING_POINTS);

		int w = n;
		for (int i = 0; i + 1 < n; i++) {
			double d = distance(path.get(i).point, path.get(i + 1).point);
			w += Math.max(0, (int) Math.ceil(d / h) - 1);
		}

		// Densified working set: current position, fixed anchor, deviation cap.
		double[] x = new double[w];
		double[] y = new double[w];
		double[] ax = new double[w];
		double[] ay = new double[w];
		double[] cap = new double[w];

		int idx = 0;
		for (int i = 0; i < n; i++) {
			TreePosition p = path.get(i);
			double clr = clearanceBound(p);
			x[idx] = p.point.x;
			y[idx] = p.point.y;
			ax[idx] = p.point.x;
			ay[idx] = p.point.y;
			cap[idx] = Math.min(tol, CLEARANCE_SAFETY * clr);
			idx++;

			if (i + 1 < n) {
				TreePosition q = path.get(i + 1);
				double clrQ = clearanceBound(q);
				double d = distance(p.point, q.point);
				int k = (int) Math.ceil(d / h) - 1;
				for (int j = 1; j <= k; j++) {
					double s = j / (double) (k + 1);
					double px = p.point.x + s * (q.point.x - p.point.x);
					double py = p.point.y + s * (q.point.y - p.point.y);
					// Lipschitz bound interpolated from both segment endpoints.
					double clrX = Math.max(0.0, Math.max(clr - s * d, clrQ - (1 - s) * d));
					x[idx] = px;
					y[idx] = py;
					ax[idx] = px;
					ay[idx] = py;
					cap[idx] = Math.min(tol, CLEARANCE_SAFETY * clrX);
					idx++;
				}
			}
		}

		double convergence = 1e-3 * h;
		for (int it = 0; it < SMOOTHING_ITERATIONS; it++) {
			double maxMove = 0;
			for (int i = 1; i < w - 1; i++) {
				double nx = x[i] + SMOOTHING_RELAX * (0.5 * (x[i - 1] + x[i + 1]) - x[i]);
				double ny = y[i] + SMOOTHING_RELAX * (0.5 * (y[i - 1] + y[i + 1]) - y[i]);

				double dx = nx - ax[i];
				double dy = ny - ay[i];
				double dd2 = dx * dx + dy * dy;
				if (dd2 > cap[i] * cap[i]) {
					double s = cap[i] / Math.sqrt(dd2);
					nx = ax[i] + dx * s;
					ny = ay[i] + dy * s;
				}

				double move = Math.abs(nx - x[i]) + Math.abs(ny - y[i]);
				if (move > maxMove) {
					maxMove = move;
				}
				x[i] = nx;
				y[i] = ny;
			}
			if (maxMove < convergence) {
				break;
			}
		}

		List<Point> out = new ArrayList<>(w);
		for (int i = 0; i < w; i++) {
			appendPoint(out, new Point(x[i], y[i]));
		}
		return out;
	}

	// -------------------------------------------------------------------------
	// Components
	// -------------------------------------------------------------------------

	private List<List<Edge>> findMedialAxisComponents() {
		List<Edge> canonical = canonicalEdges();
		Set<Edge> canonicalSet = identitySet();
		canonicalSet.addAll(canonical);

		List<List<Edge>> out = new ArrayList<>();
		Set<Edge> seen = identitySet();

		for (Edge start : canonical) {
			if (seen.contains(start)) {
				continue;
			}

			List<Edge> comp = new ArrayList<>();
			Deque<Edge> stack = new ArrayDeque<>();
			stack.push(start);

			while (!stack.isEmpty()) {
				Edge e = stack.pop();
				if (seen.contains(e)) {
					continue;
				}
				seen.add(e);
				comp.add(e);

				for (Vertex v : new Vertex[] { e.source, e.target }) {
					for (Edge oe : v.outEdges) {
						Edge c = canonicalRepresentative(oe, canonicalSet);
						if (c != null && !seen.contains(c)) {
							stack.push(c);
						}
					}
				}
			}

			out.add(comp);
		}

		return out;
	}

	private List<Edge> canonicalEdges() {
		List<Edge> out = new ArrayList<>();
		Set<Edge> seen = identitySet();

		for (Edge e : maEdges) {
			if (seen.contains(e)) {
				continue;
			}
			out.add(e);
			seen.add(e);
			if (e.twin != null) {
				seen.add(e.twin);
			}
		}
		return out;
	}

	private Edge canonicalRepresentative(Edge e, Set<Edge> canonicalSet) {
		if (e == null || !isMedialAxisEdge(e)) {
			return null;
		}
		if (canonicalSet.contains(e)) {
			return e;
		}
		if (e.twin != null && canonicalSet.contains(e.twin)) {
			return e.twin;
		}
		return null;
	}

	private boolean isMedialAxisEdge(Edge e) {
		return e.valid && e.type != EdgeType.LINESITE && e.type != EdgeType.NULLEDGE && e.type != EdgeType.OUTEDGE;
	}

	// -------------------------------------------------------------------------
	// Output conversion
	// -------------------------------------------------------------------------

	private List<List<Point>> toPointChains(List<List<TreePosition>> chains) {
		List<List<Point>> out = new ArrayList<>(chains.size());
		for (List<TreePosition> chain : chains) {
			out.add(toPoints(chain));
		}
		return out;
	}

	private List<Point> toPoints(List<TreePosition> chain) {
		List<Point> pts = new ArrayList<>(chain.size());
		for (TreePosition p : chain) {
			appendPoint(pts, p.point);
		}
		return pts;
	}

	private List<TreePosition> connectLaps(List<List<TreePosition>> laps) {
		List<TreePosition> out = new ArrayList<>();
		for (List<TreePosition> lap : laps) {
			for (TreePosition p : lap) {
				if (p.point == null) {
					continue;
				}
				if (out.isEmpty() || distance(out.get(out.size() - 1).point, p.point) > EPS) {
					out.add(p);
				}
			}
		}
		return out;
	}

	private void appendPoint(List<Point> out, Point p) {
		if (p == null) {
			return;
		}
		if (out.isEmpty() || distance(out.get(out.size() - 1), p) > EPS) {
			out.add(p);
		}
	}

	// -------------------------------------------------------------------------
	// Geometry helpers
	// -------------------------------------------------------------------------

	private void addTouchPoints(List<Point> target, Point[] pts) {
		if (pts == null) {
			return;
		}
		outer: for (Point p : pts) {
			if (p == null) {
				continue;
			}
			for (Point q : target) {
				if (distance(p, q) <= EPS) {
					continue outer;
				}
			}
			target.add(p);
		}
	}

	private Vertex otherEndpoint(Edge e, Vertex v) {
		if (e.source == v)
			return e.target;
		if (e.target == v)
			return e.source;
		return null;
	}

	private Edge orientFrom(Edge e, Vertex from) {
		if (e.source == from)
			return e;
		if (e.target == from)
			return e.twin;
		return null;
	}

	private static double distance(Point a, Point b) {
		double dx = a.x - b.x;
		double dy = a.y - b.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	private static Point interpolate(Point a, Point b, double t) {
		return new Point(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
	}

	private static double clamp(double x, double lo, double hi) {
		return Math.max(lo, Math.min(hi, x));
	}

	private static double angleOf(Point from, Point to) {
		return FastMath.atan2(to.y - from.y, to.x - from.x);
	}

	private static double ccwDelta(double base, double angle) {
		double d = (angle - base) % (2 * Math.PI);
		return d < 0 ? d + 2 * Math.PI : d;
	}

	private static <T> Set<T> identitySet() {
		return Collections.newSetFromMap(new IdentityHashMap<>());
	}

	// -------------------------------------------------------------------------
	// Internal data
	// -------------------------------------------------------------------------

	private static class DiscreteTree {
		DiscreteNode root;
		final List<DiscreteNode> maNodes = new ArrayList<>();
	}

	private static class DiscreteNode {
		Point point;
		double clearance;
		boolean boundaryLeaf;

		final List<Point> touchPoints = new ArrayList<>();

		DiscreteNode parent;
		DiscreteEdge parentEdge;
		final List<DiscreteEdge> children = new ArrayList<>();

		double height;
		double startTime;

		DiscreteEdge bestChild;
		DiscreteNode bestLeaf;

		// binary-lifting ancestor table: up[k] is the 2^k-th ancestor
		int depth;
		DiscreteNode[] up;

		Branch branch; // only used on leaves
	}

	private static class DiscreteEdge {
		DiscreteNode parent;
		DiscreteNode child;
		double length;
		double sortKey; // scratch for sortChildrenByEmbedding
	}

	private static class Branch {
		int index;
		DiscreteNode root;
		DiscreteNode leaf;
	}

	private static class TreePosition {
		Branch branch;
		Point point;
		double time;
		double height;

		DiscreteNode node;
		DiscreteEdge edge;
		double offset;
	}
}
