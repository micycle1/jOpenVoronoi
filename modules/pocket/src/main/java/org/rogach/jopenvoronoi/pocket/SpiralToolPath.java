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

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Raw spiral tool-path generation based on the discrete medial-axis tree
 * approach of Held / Spielberger.
 *
 * <p>
 * This implementation is intentionally conservative:
 * <ul>
 * <li>it roots the medial axis at the largest-clearance vertex,</li>
 * <li>builds a rooted tree of medial-axis edges,</li>
 * <li>discretizes each non-linear MA edge into a fixed number of samples,</li>
 * <li>adds clearance-line leaves from {@link Edge#micTouchPoints(double)},</li>
 * <li>computes impulse start-times / heights,</li>
 * <li>builds wavefronts and laps,</li>
 * <li>connects laps into a raw spiral polyline.</li>
 * </ul>
 *
 * <p>
 * Intended for simple pockets (tree-like medial axis). Islands / cyclic MA
 * components are not supported here.
 */
public class SpiralToolPath {

	private static final double EPS = 1e-9;

	private final HalfEdgeDiagram g;
	private final List<Edge> maEdges = new ArrayList<>();

	private double maxStepOver = 0.05;
	private int samplesPerCurveEdge = 10;
	private int maxLaps = 2000;

	private boolean storeWavefronts = false;
	private boolean storeLaps = false;

	private final List<List<Point>> toolPathComponents = new ArrayList<>();
	private final List<List<List<Point>>> lapComponents = new ArrayList<>();
	private final List<List<List<Point>>> wavefrontComponents = new ArrayList<>();

	public SpiralToolPath(HalfEdgeDiagram g) {
		this.g = g;
		for (Edge e : g.edges) {
			if (isMedialAxisEdge(e)) {
				maEdges.add(e);
			}
		}
	}

	public void setStepOver(double stepOver) {
		if (stepOver <= 0) {
			throw new IllegalArgumentException("stepOver must be > 0");
		}
		this.maxStepOver = stepOver;
	}

	public void setWidth(double stepOver) {
		setStepOver(stepOver);
	}

	/**
	 * Number of uniform samples per non-linear MA edge. A value of 10 means 10
	 * segments / 11 sampled points including endpoints. Linear edges are not
	 * subdivided.
	 */
	public void setSamplesPerCurveEdge(int n) {
		if (n < 1) {
			throw new IllegalArgumentException("samplesPerCurveEdge must be >= 1");
		}
		this.samplesPerCurveEdge = n;
	}

	public void setMaxLaps(int maxLaps) {
		if (maxLaps < 1) {
			throw new IllegalArgumentException("maxLaps must be >= 1");
		}
		this.maxLaps = maxLaps;
	}

	public void setStoreWavefronts(boolean storeWavefronts) {
		this.storeWavefronts = storeWavefronts;
	}

	public void setStoreLaps(boolean storeLaps) {
		this.storeLaps = storeLaps;
	}

	/**
	 * Returns the generated raw spiral tool path for each connected medial-axis
	 * component.
	 * <p>
	 * Each returned component is a single polyline obtained by concatenating the
	 * successive laps produced by the spiral algorithm, starting near the root
	 * (largest-clearance point) and progressing outward toward the pocket boundary.
	 * <p>
	 * This data is always populated after a successful call to {@link #run()}.
	 *
	 * @return a list of raw spiral polylines, one per connected medial-axis
	 *         component
	 */
	public List<List<Point>> getToolPathComponents() {
		return toolPathComponents;
	}

	/**
	 * Returns the lap geometry generated for each connected medial-axis component.
	 * <p>
	 * A lap corresponds to one closed-like outward pass of the raw spiral
	 * construction. The returned nesting is:
	 * 
	 * <pre>
	 * component -&gt; lap -&gt; polyline points
	 * </pre>
	 * 
	 * Lap {@code 0} is the innermost lap and the last lap lies nearest the pocket
	 * boundary.
	 * <p>
	 * Lap data is only retained when enabled via {@link #setStoreLaps(boolean)}. If
	 * storage is disabled, this method returns an empty list even after
	 * {@link #run()}.
	 *
	 * @return lap polylines per component, or an empty list when lap storage is
	 *         disabled
	 */
	public List<List<List<Point>>> getLapComponents() {
		return lapComponents;
	}

	/**
	 * Returns the sampled wavefronts generated for each connected medial-axis
	 * component.
	 * <p>
	 * A wavefront is the polygonal chain obtained by evaluating the discrete
	 * impulse-propagation model at a specific time value. The returned nesting is:
	 * 
	 * <pre>
	 * component -&gt; wavefront -&gt; polyline points
	 * </pre>
	 * 
	 * Wavefronts are ordered by increasing time, from the initial wavefront near
	 * the root to the final wavefront at the pocket boundary.
	 * <p>
	 * Wavefront data is only retained when enabled via
	 * {@link #setStoreWavefronts(boolean)}. If storage is disabled, this method
	 * returns an empty list even after {@link #run()}.
	 *
	 * @return wavefront polylines per component, or an empty list when wavefront
	 *         storage is disabled
	 */
	public List<List<List<Point>>> getWavefrontComponents() {
		return wavefrontComponents;
	}

	public void run() {
		toolPathComponents.clear();
		lapComponents.clear();
		wavefrontComponents.clear();

		List<List<Edge>> components = findMedialAxisComponents();
		for (List<Edge> component : components) {
			if (component.isEmpty()) {
				continue;
			}

			Set<Edge> componentCanonicalSet = identitySet();
			componentCanonicalSet.addAll(component);

			DiscreteTree discreteTree = buildDiscreteTree(component, componentCanonicalSet);
			if (discreteTree == null || discreteTree.root == null) {
				continue;
			}
			if (discreteTree == null || discreteTree.root == null) {
				continue;
			}

			sortChildrenByEmbedding(discreteTree.root);
			computeHeightsAndStartTimes(discreteTree.root);

			List<DiscreteNode> leafOrder = new ArrayList<>();
			collectLeaves(discreteTree.root, leafOrder);
			if (leafOrder.isEmpty()) {
				continue;
			}

			List<Branch> branches = buildBranches(leafOrder);

			int m = Math.max(1, (int) Math.ceil(discreteTree.root.height / maxStepOver));
			if (m > maxLaps) {
				throw new IllegalArgumentException(
						"Requested step-over creates too many laps: " + m + " (height=" + discreteTree.root.height + ", stepOver=" + maxStepOver + ")");
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

			toolPathComponents.add(connectLaps(laps));
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

	private DiscreteTree buildDiscreteTree(List<Edge> component, Set<Edge> componentCanonicalSet) {
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
				int segments = segmentCount(oriented);

				for (int i = 1; i < segments; i++) {
					double u = i / (double) segments;
					DiscreteNode sampleNode = createSampleNode(oriented, u, tree);
					prev = connectOrMerge(prev, sampleNode, tree);
				}

				DiscreteNode childVertexNode = createVertexNode(neighbor, componentCanonicalSet, tree);
				prev = connectOrMerge(prev, childVertexNode, tree);

				visited.add(neighbor);
				stack.push(new VertexFrame(neighbor, childVertexNode));
			}
		}

		attachClearanceLeaves(tree);
		return tree;
	}

	private int segmentCount(Edge e) {
		return (e.type == EdgeType.LINE) ? 1 : samplesPerCurveEdge;
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
		n.id = tree.nextId++;
		n.point = v.position;
		n.clearance = v.dist();

		for (Edge oe : v.outEdges) {
			Edge canonical = canonicalRepresentative(oe, componentCanonicalSet);
			if (canonical == null) {
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
		n.id = tree.nextId++;
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
			tree.maNodes.remove(child);
			return parent;
		}

		DiscreteEdge de = new DiscreteEdge();
		de.parent = parent;
		de.child = child;
		de.length = len;

		child.parent = parent;
		child.parentEdge = de;
		parent.children.add(de);
		tree.edges.add(de);

		return child;
	}

	private void attachClearanceLeaves(DiscreteTree tree) {
		for (DiscreteNode center : new ArrayList<>(tree.maNodes)) {
			for (Point tp : center.touchPoints) {
				if (tp == null || distance(center.point, tp) <= EPS) {
					continue;
				}

				DiscreteNode leaf = new DiscreteNode();
				leaf.id = tree.nextId++;
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
				tree.edges.add(edge);
			}
		}
	}

	// -------------------------------------------------------------------------
	// Tree ordering / height / times
	// -------------------------------------------------------------------------

	private void sortChildrenByEmbedding(DiscreteNode node) {
		if (node.children.size() > 1) {
			if (node.parent == null) {
				node.children.sort(Comparator.comparingDouble(e -> angleOf(node.point, e.child.point)));
			} else {
				final double base = angleOf(node.point, node.parent.point);
				node.children.sort(Comparator.comparingDouble(e -> ccwDelta(base, angleOf(node.point, e.child.point))));
			}
		}

		for (DiscreteEdge e : node.children) {
			sortChildrenByEmbedding(e.child);
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

	private void collectLeaves(DiscreteNode node, List<DiscreteNode> out) {
		if (node.children.isEmpty()) {
			out.add(node);
			return;
		}
		for (DiscreteEdge e : node.children) {
			collectLeaves(e.child, out);
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

	private List<Branch> buildBranches(List<DiscreteNode> leafOrder) {
		List<Branch> out = new ArrayList<>(leafOrder.size());

		for (int i = 0; i < leafOrder.size(); i++) {
			DiscreteNode leaf = leafOrder.get(i);

			List<DiscreteEdge> revEdges = new ArrayList<>();
			DiscreteNode root = leaf;

			while (root.parent != null) {
				revEdges.add(root.parentEdge);
				root = root.parent;
			}
			Collections.reverse(revEdges);

			Branch b = new Branch();
			b.index = i;
			b.root = root;
			b.leaf = leaf;
			b.edges = revEdges.toArray(new DiscreteEdge[0]);
			b.childStartTimes = new double[b.edges.length];

			for (int k = 0; k < b.edges.length; k++) {
				b.childStartTimes[k] = b.edges[k].child.startTime;
			}

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

		if (branch.edges.length == 0) {
			return positionAtNode(branch, branch.leaf);
		}

		if (time <= branch.root.startTime + EPS) {
			return positionAtNode(branch, branch.root);
		}
		if (time >= branch.leaf.startTime - EPS) {
			return positionAtNode(branch, branch.leaf);
		}

		int idx = firstEdgeEndingAfter(branch.childStartTimes, time);
		if (idx >= branch.edges.length) {
			return positionAtNode(branch, branch.leaf);
		}

		DiscreteEdge e = branch.edges[idx];
		double t0 = e.parent.startTime;
		double t1 = e.child.startTime;
		double alpha = (t1 - t0 <= EPS) ? 1.0 : clamp((time - t0) / (t1 - t0), 0.0, 1.0);

		return positionOnEdge(branch, e, alpha * e.length, time);
	}

	private int firstEdgeEndingAfter(double[] times, double time) {
		int lo = 0, hi = times.length;
		while (lo < hi) {
			int mid = (lo + hi) >>> 1;
			if (times[mid] <= time + EPS) {
				lo = mid + 1;
			} else {
				hi = mid;
			}
		}
		return lo;
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
			List<Point> pts = new ArrayList<>(chain.size());
			for (TreePosition p : chain) {
				appendPoint(pts, p.point);
			}
			out.add(pts);
		}
		return out;
	}

	private List<Point> connectLaps(List<List<TreePosition>> laps) {
		List<Point> out = new ArrayList<>();
		for (List<TreePosition> lap : laps) {
			for (TreePosition p : lap) {
				appendPoint(out, p.point);
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
		return a.sub(b).norm();
	}

	private static Point interpolate(Point a, Point b, double t) {
		return a.add(b.sub(a).mult(t));
	}

	private static double clamp(double x, double lo, double hi) {
		return Math.max(lo, Math.min(hi, x));
	}

	private static double angleOf(Point from, Point to) {
		return Math.atan2(to.y - from.y, to.x - from.x);
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
		int nextId = 0;
		DiscreteNode root;
		final List<DiscreteNode> maNodes = new ArrayList<>();
		final List<DiscreteEdge> edges = new ArrayList<>();
	}

	private static class DiscreteNode {
		int id;
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

		Branch branch; // only used on leaves
	}

	private static class DiscreteEdge {
		DiscreteNode parent;
		DiscreteNode child;
		double length;
	}

	private static class Branch {
		int index;
		DiscreteNode root;
		DiscreteNode leaf;
		DiscreteEdge[] edges;
		double[] childStartTimes;
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