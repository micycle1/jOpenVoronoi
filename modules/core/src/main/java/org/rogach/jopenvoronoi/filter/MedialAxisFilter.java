package org.rogach.jopenvoronoi.filter;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;
import org.rogach.jopenvoronoi.vertex.VertexType;

/**
 * Prunes low-salience medial-axis branches that are dominated by near-parallel
 * generating segments (e.g. the spurious spikes a medial axis grows opposite a
 * sharp reflex corner), while always keeping site/null edges and cyclic medial
 * branches.
 * <p>
 * All branch keep/discard decisions are computed once, eagerly, in
 * {@link #setGraph(HalfEdgeDiagram)} — before {@link org.rogach.jopenvoronoi.VoronoiDiagram#filter}
 * starts clearing {@code valid} flags on this same pass. This makes the result
 * independent of {@code g.edges} iteration order: {@link #apply(Edge)} is a
 * pure lookup into a decision table snapshotted from the diagram's state at the
 * start of the pass, rather than a lazy, order-sensitive graph walk. Recomputing
 * the table on every {@link #setGraph} call also makes filter instances safe to
 * reuse across diagrams, or across repeated filtering passes on a diagram that
 * has changed since the last pass.
 */
public class MedialAxisFilter extends Filter {

	/**
	 * Dot-product threshold in [0,1] for treating two line-site segments as nearly
	 * parallel. Comparison uses abs(dot), so opposite directions still count as
	 * parallel.
	 */
	private final double dotProductThreshold;

	/**
	 * Fraction of branch edges that must look near-parallel before the branch is
	 * considered "parallel-dominated".
	 */
	private final double parallelFractionThreshold;

	/**
	 * A branch is considered length-salient if
	 *
	 * branchLength >= minLengthFactor * maxRadius
	 */
	private final double minLengthFactor;

	/**
	 * A branch is considered radius-salient if
	 *
	 * radiusGain >= minRadiusGainFactor * maxRadius
	 *
	 * where radiusGain = maxRadiusAlongBranch - attachmentRadius.
	 */
	private final double minRadiusGainFactor;

	private static final double EPS = 1e-9;

	/**
	 * Keep/discard decision for every edge in the diagram, computed fresh by
	 * {@link #setGraph(HalfEdgeDiagram)} on every filtering pass.
	 */
	private Map<Edge, Boolean> decisions = Collections.emptyMap();

	public MedialAxisFilter() {
		this(0.8, 0.5, 2.0, 0.25);
	}

	/**
	 * Creates a medial-axis filter with a custom near-parallel threshold and
	 * default branch-pruning parameters.
	 * <p>
	 * The threshold is applied to the absolute dot product of the two adjacent
	 * line-site directions. Values closer to {@code 1} treat only more strongly
	 * aligned segments as parallel, while lower values classify a wider range of
	 * segment pairs as parallel.
	 * <p>
	 * This constructor uses the following defaults:
	 * <ul>
	 * <li>{@code parallelFractionThreshold = 0.5}</li>
	 * <li>{@code minLengthFactor = 2.0}</li>
	 * <li>{@code minRadiusGainFactor = 0.25}</li>
	 * </ul>
	 *
	 * @param dotProductThreshold threshold in {@code [0,1]} for treating two
	 *                            adjacent line-site segments as nearly parallel,
	 *                            based on {@code abs(dot)}
	 */
	public MedialAxisFilter(double dotProductThreshold) {
		this(dotProductThreshold, 0.5, 2.0, 0.25);
	}

	/**
	 * Creates a medial-axis filter with fully configurable branch-pruning
	 * parameters.
	 * <p>
	 * A medial branch is pruned only when it is both:
	 * <ul>
	 * <li>dominated by near-parallel generating segments, and</li>
	 * <li>insignificant according to branch salience tests</li>
	 * </ul>
	 * <p>
	 * The parameters are interpreted as follows:
	 * <ul>
	 * <li>{@code dotProductThreshold}: threshold in {@code [0,1]} applied to
	 * {@code abs(dot)} of the two line-site directions; larger values prune less
	 * aggressively</li>
	 * <li>{@code parallelFractionThreshold}: minimum fraction of tested edges in a
	 * branch that must be classified as near-parallel before the whole branch is
	 * considered parallel-dominated</li>
	 * <li>{@code minLengthFactor}: a branch is considered salient by length when
	 * {@code branchLength >= minLengthFactor * maxRadius}</li>
	 * <li>{@code minRadiusGainFactor}: a branch is considered salient by radius
	 * when {@code radiusGain >= minRadiusGainFactor * maxRadius}, where
	 * {@code radiusGain = maxRadiusAlongBranch - attachmentRadius}</li>
	 * </ul>
	 *
	 * @param dotProductThreshold       threshold in {@code [0,1]} for classifying
	 *                                  two adjacent line-site segments as nearly
	 *                                  parallel
	 * @param parallelFractionThreshold minimum fraction in {@code [0,1]} of branch
	 *                                  edges that must test as near-parallel before
	 *                                  the branch is treated as parallel-dominated
	 * @param minLengthFactor           minimum normalized branch-length factor
	 *                                  required for length-based salience
	 * @param minRadiusGainFactor       minimum normalized radius-gain factor
	 *                                  required for radius-based salience
	 */
	public MedialAxisFilter(double dotProductThreshold, double parallelFractionThreshold, double minLengthFactor, double minRadiusGainFactor) {
		this.dotProductThreshold = dotProductThreshold;
		this.parallelFractionThreshold = parallelFractionThreshold;
		this.minLengthFactor = minLengthFactor;
		this.minRadiusGainFactor = minRadiusGainFactor;
	}

	/**
	 * Snapshots the diagram and eagerly computes a keep/discard decision for every
	 * edge, before any edge in this pass has had its {@code valid} flag mutated.
	 */
	@Override
	public void setGraph(HalfEdgeDiagram g) {
		super.setGraph(g);
		decisions = computeDecisions();
	}

	/**
	 * Looks up the decision computed for {@code e} by {@link #setGraph}. Never
	 * walks the graph itself, so the result cannot depend on the order in which
	 * {@link org.rogach.jopenvoronoi.VoronoiDiagram#filter} visits edges.
	 */
	@Override
	public boolean apply(Edge e) {
		if (e == null) {
			return false;
		}
		Boolean decision = decisions.get(e);
		return decision != null && decision.booleanValue();
	}

	// -------------------------------------------------------------------------
	// Eager decision computation
	// -------------------------------------------------------------------------

	private Map<Edge, Boolean> computeDecisions() {
		Map<Edge, Boolean> result = new IdentityHashMap<>();
		if (g == null) {
			return result;
		}
		Set<Edge> assigned = Collections.newSetFromMap(new IdentityHashMap<>());

		for (Edge e : g.edges) {
			if (result.containsKey(e)) {
				continue;
			}
			if (!e.valid) {
				result.put(e, Boolean.FALSE);
				continue;
			}
			if (e.type == EdgeType.LINESITE || e.type == EdgeType.NULLEDGE) {
				result.put(e, Boolean.TRUE);
				continue;
			}
			if (e.type == EdgeType.SEPARATOR) {
				result.put(e, Boolean.FALSE);
				continue;
			}
			if (assigned.contains(e)) {
				continue; // already covered by a branch built from a previously visited edge
			}
			BranchInfo b = buildBranch(e, assigned);
			for (Edge be : b.edges) {
				result.put(be, b.keep);
				if (be.twin != null) {
					result.put(be.twin, b.keep);
				}
			}
		}
		return result;
	}

	// -------------------------------------------------------------------------
	// Branch extraction / evaluation
	// -------------------------------------------------------------------------

	private BranchInfo buildBranch(Edge seed, Set<Edge> assigned) {
		Deque<Edge> branchEdges = new ArrayDeque<>();

		// Store one outgoing half-edge per undirected edge, ordered along branch.
		branchEdges.add(seed);
		markAssigned(seed, assigned);

		Vertex startVertex = extendBackward(branchEdges, seed.source, seed, assigned);
		Vertex endVertex = extendForward(branchEdges, seed.target, seed.twin, assigned);

		boolean cycle = (medialDegree(startVertex) == 2) && (medialDegree(endVertex) == 2);

		List<Edge> orderedEdges = new ArrayList<>(branchEdges);
		return evaluateBranch(orderedEdges, startVertex, endVertex, cycle);
	}

	/**
	 * Walk from the seed's source side until a non-degree-2 vertex or a cycle
	 * closure is reached.
	 *
	 * @param currentVertex   current branch vertex
	 * @param currentAtVertex outgoing half-edge from currentVertex representing the
	 *                        current undirected edge
	 */
	private Vertex extendBackward(Deque<Edge> branchEdges, Vertex currentVertex, Edge currentAtVertex, Set<Edge> assigned) {
		while (medialDegree(currentVertex) == 2) {
			Edge next = otherMedialEdge(currentVertex, currentAtVertex);
			if (next == null || assigned.contains(next)) {
				break;
			}

			branchEdges.addFirst(next);
			markAssigned(next, assigned);

			currentVertex = next.target;
			currentAtVertex = next.twin;
		}
		return currentVertex;
	}

	/**
	 * Walk from the seed's target side until a non-degree-2 vertex or a cycle
	 * closure is reached.
	 *
	 * @param currentVertex   current branch vertex
	 * @param currentAtVertex outgoing half-edge from currentVertex representing the
	 *                        current undirected edge
	 */
	private Vertex extendForward(Deque<Edge> branchEdges, Vertex currentVertex, Edge currentAtVertex, Set<Edge> assigned) {
		while (medialDegree(currentVertex) == 2) {
			Edge next = otherMedialEdge(currentVertex, currentAtVertex);
			if (next == null || assigned.contains(next)) {
				break;
			}

			branchEdges.addLast(next);
			markAssigned(next, assigned);

			currentVertex = next.target;
			currentAtVertex = next.twin;
		}
		return currentVertex;
	}

	private void markAssigned(Edge e, Set<Edge> assigned) {
		assigned.add(e);
		if (e.twin != null) {
			assigned.add(e.twin);
		}
	}

	private int medialDegree(Vertex v) {
		int degree = 0;
		for (Edge e : v.outEdges) {
			if (isMedialEdge(e)) {
				degree++;
			}
		}
		return degree;
	}

	/**
	 * Returns the other outgoing medial half-edge at the vertex, excluding the
	 * outgoing half-edge that represents the current undirected edge there.
	 */
	private Edge otherMedialEdge(Vertex v, Edge exclude) {
		for (Edge e : v.outEdges) {
			if (e == exclude) {
				continue;
			}
			if (isMedialEdge(e)) {
				return e;
			}
		}
		return null;
	}

	private boolean isMedialEdge(Edge e) {
		return e != null && e.valid && e.type != EdgeType.LINESITE && e.type != EdgeType.NULLEDGE && e.type != EdgeType.SEPARATOR;
	}

	private BranchInfo evaluateBranch(List<Edge> edges, Vertex startVertex, Vertex endVertex, boolean cycle) {
		double branchLength = 0.0;
		double maxRadius = 0.0;

		for (Edge e : edges) {
			branchLength += safeLength(e);
			maxRadius = Math.max(maxRadius, radius(e.source));
			maxRadius = Math.max(maxRadius, radius(e.target));
		}

		double startRadius = radius(startVertex);
		double endRadius = radius(endVertex);

		// For a leaf-junction branch, attachment radius should be the trunk side,
		// which is typically the larger endpoint radius.
		double attachmentRadius = cycle ? maxRadius : Math.max(startRadius, endRadius);
		double radiusGain = Math.max(0.0, maxRadius - attachmentRadius);

		int parallelCount = 0;
		int testedCount = 0;

		for (Edge e : edges) {
			Boolean nearParallel = segmentsParallelSafe(e);
			if (nearParallel != null) {
				testedCount++;
				if (nearParallel.booleanValue()) {
					parallelCount++;
				}
			}
		}

		double parallelFraction = testedCount == 0 ? 0.0 : ((double) parallelCount) / testedCount;

		boolean salientByLength = branchLength >= minLengthFactor * Math.max(maxRadius, EPS);

		boolean salientByRadius = radiusGain >= minRadiusGainFactor * Math.max(maxRadius, EPS);

		boolean parallelDominated = parallelFraction >= parallelFractionThreshold;

		// near-parallel is only a penalty, not an automatic delete.
		// We prune only if the branch is low-salience AND parallel-dominated.
		boolean keep = cycle || salientByLength || salientByRadius || !parallelDominated;

		return new BranchInfo(edges, keep);
	}

	private double safeLength(Edge e) {
		double len = e.length();
		return Double.isFinite(len) && len > 0.0 ? len : 0.0;
	}

	private double radius(Vertex v) {
		if (v == null) {
			return 0.0;
		}
		double d = v.dist();
		return Double.isFinite(d) && d > 0.0 ? d : 0.0;
	}

	// -------------------------------------------------------------------------
	// Local near-parallel test
	// -------------------------------------------------------------------------

	/**
	 * @return TRUE/FALSE if the test was possible, or {@code null} if the local
	 *         topology did not match the expected endpoint pattern. Unlike a
	 *         broad catch-all, an unexpected failure elsewhere (e.g. a bug
	 *         surfacing as an NPE) is not swallowed here and propagates normally.
	 */
	private Boolean segmentsParallelSafe(Edge e) {
		Vertex endp1 = findEndpoint(e);
		Vertex endp2 = findEndpoint(e.twin);
		if (endp1 == null || endp2 == null) {
			return null;
		}

		Edge e1 = findSegment(endp1);
		Edge e2 = findSegment(endp2);
		if (e1 == null || e2 == null) {
			return null;
		}
		e2 = e2.twin;

		double dotprod = Math.abs(edgeDotprod(e1, e2));
		return Boolean.valueOf(dotprod > dotProductThreshold);
	}

	/**
	 * Calculate the dot-product between unit vectors aligned along edges e1 -> e2.
	 */
	private double edgeDotprod(Edge e1, Edge e2) {
		// guard against zero-length site segments before normalization.
		if (safeLength(e1) <= EPS || safeLength(e2) <= EPS) {
			return 0.0;
		}

		Vertex src1 = e1.source;
		Vertex trg1 = e1.target;
		Vertex src2 = e2.source;
		Vertex trg2 = e2.target;

		Point sp1 = src1.position;
		Point tp1 = trg1.position;
		Point sp2 = src2.position;
		Point tp2 = trg2.position;

		Point dir1 = tp1.sub(sp1);
		Point dir2 = tp2.sub(sp2);

		dir1.normalize();
		dir2.normalize();

		return dir1.dot(dir2);
	}

	/**
	 * Finds the line-site edge incident to the given endpoint vertex.
	 *
	 * @param v endpoint vertex on the polygon boundary
	 * @return incident line-site edge, or {@code null} if none is found
	 */
	private Edge findSegment(Vertex v) {
		if (v == null) {
			return null;
		}
		for (Edge e : v.outEdges) {
			if (e.type == EdgeType.LINESITE) {
				return e;
			}
		}
		return null;
	}

	/**
	 * Finds the endpoint vertex reached through a neighboring null edge.
	 *
	 * @param e edge from the Voronoi diagram
	 * @return polygon endpoint associated with e, or {@code null} if the
	 *         surrounding topology doesn't match the expected pattern
	 */
	private Vertex findEndpoint(Edge e) {
		if (e == null) {
			return null;
		}
		Edge next = e.next;
		Edge prev = e.prev;

		if (next != null && next.type == EdgeType.NULLEDGE) {
			Vertex endp = next.target;
			assert (endp.type == VertexType.ENDPOINT) : "endp.type == VertexType.ENDPOINT";
			return endp;
		}

		if (prev != null && prev.type == EdgeType.NULLEDGE) {
			Vertex endp = prev.source;
			assert (endp.type == VertexType.ENDPOINT) : "endp.type == VertexType.ENDPOINT";
			return endp;
		}

		return null;
	}

	// -------------------------------------------------------------------------
	// Branch data
	// -------------------------------------------------------------------------

	private static final class BranchInfo {
		final List<Edge> edges;
		final boolean keep;

		BranchInfo(List<Edge> edges, boolean keep) {
			this.edges = edges;
			this.keep = keep;
		}
	}
}
