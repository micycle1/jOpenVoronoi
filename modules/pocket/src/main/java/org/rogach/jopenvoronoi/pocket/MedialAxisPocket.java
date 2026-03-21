package org.rogach.jopenvoronoi.pocket;

import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.ArrayDeque;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AllowedSolution;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Medial-axis pocketing: walks the medial axis of a Voronoi diagram and
 * produces a sequence of {@link MIC Maximal Inscribed Circles} whose
 * envelopes cover the pocket interior.
 *
 * <p><b>Algorithm overview</b></p>
 * <ol>
 *   <li>The input graph must be a Voronoi diagram filtered through
 *       {@link org.rogach.jopenvoronoi.filter.PolygonInteriorFilter PolygonInteriorFilter}
 *       and then {@link org.rogach.jopenvoronoi.filter.MedialAxisFilter MedialAxisFilter}
 *       so that only interior medial-axis edges remain.</li>
 *   <li>The algorithm begins at the vertex with the <em>largest</em>
 *       clearance-disk radius (the deepest interior point) and walks outward
 *       along medial-axis edges, emitting MICs at intervals governed by
 *       {@link #setWidth(double) maxWidth}.</li>
 *   <li>At medial-axis junctions the algorithm pushes un-visited branches
 *       onto a stack and backtracks to them once the current branch ends.</li>
 *   <li>Separate connected components of the medial axis are processed
 *       independently.</li>
 * </ol>
 *
 * <p><b>Cut width</b></p>
 * <p>
 * The <em>cut width</em> (also called <em>step-over</em> in CNC terminology)
 * is the maximum width of new material removed by a single MIC step. It is
 * defined as:
 * <pre>
 *   w = |c₂ − c₁| + r₂ − r₁
 * </pre>
 * where (c₁, r₁) is the previously-cut circle and (c₂, r₂) is the new one.
 * <p>
 * A smaller cut width produces more, closely-spaced MICs (finer step-over,
 * better surface finish, lower cutting forces). A larger cut width produces
 * fewer MICs but may still generate small circles in convex corners where the
 * medial axis tapers to a point — the algorithm simply reaches those corners
 * in fewer steps.
 *
 * <p><b>Connected components</b></p>
 * <p>
 * {@link #getMicComponents()} returns one MIC list per connected component of
 * the medial axis. A simple polygon (no holes) produces a single connected
 * medial axis and therefore one component. When the pocket contains internal
 * holes (islands), the medial axis splits into multiple disconnected
 * components, each of which is processed independently and returned as a
 * separate list.
 *
 * <p><b>Usage</b></p>
 * <pre>
 * {@code
 * VoronoiDiagram vd = new VoronoiDiagram();
 * // ... insert point/line sites ...
 * vd.filter(new PolygonInteriorFilter(true));
 * vd.filter(new MedialAxisFilter());
 *
 * MedialAxisPocket pocket = new MedialAxisPocket(vd.getDiagram());
 * pocket.setWidth(0.05);   // 0.05 units step-over
 * pocket.run();
 * List<List<MIC>> components = pocket.getMicComponents();
 * }
 * </pre>
 *
 * @see MIC
 */
public class MedialAxisPocket {

	private final HalfEdgeDiagram g;
	private final List<Edge> maEdges = new ArrayList<>();
	private final Map<Edge, Boolean> edgeDone = new HashMap<>();
	private final Deque<BranchPoint> unvisited = new ArrayDeque<>();
	private Edge currentEdge;
	private double currentRadius;
	private double currentU;
	private Point currentCenter;
	private boolean newBranch;
	private Point previousBranchCenter;
	private double previousBranchRadius;
	private double maxWidth;
	private List<MIC> micList = new ArrayList<>();
	private final List<List<MIC>> maComponents = new ArrayList<>();

	/**
	 * Creates a new medial-axis pocket processor.
	 *
	 * @param g a Voronoi diagram graph that has been filtered through a
	 *          {@link org.rogach.jopenvoronoi.filter.MedialAxisFilter MedialAxisFilter}
	 */
	public MedialAxisPocket(HalfEdgeDiagram g) {
		this.g = g;
		for (Edge e : g.edges) {
			if (e.valid
					&& e.type != EdgeType.LINESITE
					&& e.type != EdgeType.NULLEDGE
					&& e.type != EdgeType.OUTEDGE) {
				maEdges.add(e);
				edgeDone.put(e, false);
			}
		}
		currentEdge = null;
		maxWidth = 0.05;
	}

	/**
	 * Sets the maximum cut width (step-over).
	 * <p>
	 * In CNC milling this controls how much fresh material is engaged on each
	 * step: a smaller value yields a finer finish but more MICs (and therefore
	 * more tool moves), while a larger value is more aggressive but still
	 * produces small MICs where the medial axis tapers into convex corners.
	 *
	 * @param w maximum cut width in diagram units (must be positive)
	 */
	public void setWidth(double w) {
		maxWidth = w;
	}

	/**
	 * Returns the algorithm output: one MIC list per connected component of the
	 * medial axis.
	 * <p>
	 * A simple polygon (no holes) has a single connected medial axis and will
	 * produce exactly one component. Pockets with interior holes (islands) cause
	 * the medial axis to split into multiple disconnected components, each
	 * returned as a separate list.
	 *
	 * @return list of MIC lists, one per connected component
	 */
	public List<List<MIC>> getMicComponents() {
		return maComponents;
	}

	/**
	 * Run the pocketing algorithm across all connected components of the medial
	 * axis.
	 */
	public void run() {
		micList.clear();
		while (findInitialMic()) {
			while (findNextMic()) {
				// keep finding MICs on this component
			}
			maComponents.add(micList);
			micList = new ArrayList<>();
		}
	}

	/**
	 * Find the largest MIC across all undone medial-axis edges and set it as
	 * the starting point. Returns {@code false} when no undone edges remain.
	 */
	private boolean findInitialMic() {
		MIC mic = new MIC();
		double maxMicRadius = -1;
		Point maxMicPos = new Point(0, 0);
		Vertex maxMicVertex = null;
		boolean found = false;

		for (Edge e : maEdges) {
			Vertex src = e.source;
			if (!edgeDone.get(e) && src.dist() > maxMicRadius) {
				maxMicRadius = src.dist();
				maxMicPos = src.position;
				maxMicVertex = src;
				found = true;
			}
		}
		if (!found) {
			return false;
		}

		mic.c2 = maxMicPos;
		mic.r2 = maxMicRadius;
		mic.c1 = mic.c2;
		mic.r1 = mic.r2;
		currentU = 0;
		currentRadius = maxMicRadius;
		currentCenter = maxMicPos;
		previousBranchCenter = maxMicPos;
		previousBranchRadius = maxMicRadius;
		mic.cPrev = maxMicPos;
		mic.rPrev = maxMicRadius;

		// find the edge on which we start machining: the out-edge whose target
		// has the largest clearance-disk radius
		double maxAdjRadius = -1;
		for (Edge e : maxMicVertex.outEdges) {
			if (e.target.dist() > maxAdjRadius && e.valid && e.type != EdgeType.OUTEDGE) {
				maxAdjRadius = e.target.dist();
				currentEdge = e;
			}
		}
		// stash the other valid medial-axis out-edges for visiting later
		for (Edge e : maxMicVertex.outEdges) {
			if (e != currentEdge && e.valid
					&& e.type != EdgeType.OUTEDGE
					&& e.type != EdgeType.NULLEDGE
					&& e.type != EdgeType.LINESITE) {
				unvisited.push(new BranchPoint(currentCenter, currentRadius, e));
			}
		}

		newBranch = false;
		mic.newBranch = newBranch;
		micList.add(mic);
		
		Point[] fp = edgeFootPoints(currentEdge, 0.0);
		if (fp != null && fp.length == 2) {
			mic.tp1a = fp[0];
			mic.tp1b = fp[1];
			mic.tp2a = fp[0];
			mic.tp2b = fp[1];
		}
		return true;
	}

	/**
	 * Attempt to find and output the next MIC. Returns {@code true} if a MIC
	 * was added, {@code false} for end-of-operation.
	 */
	private boolean findNextMic() {
		if (currentEdge == null) {
			return false;
		}

		Entry<Point, Double> targetSample = edgePoint(currentEdge, 1.0);
		Point c2 = targetSample.getKey();
		double r2 = targetSample.getValue();
		double wTarget = cutWidth(currentCenter, currentRadius, c2, r2);

		if (wTarget > maxWidth) {
			// moving to the target vertex would exceed cut-width →
			// search along the current edge for the next MIC
			Entry<Double, Double> nextUR = findNextU();
			outputNextMic(nextUR.getKey(), nextUR.getValue(), newBranch);
			return true;
		} else {
			// moving to the target is within cut-width → advance to a new edge
			markDone(currentEdge);

			Entry<Edge, Boolean> nextEdgeResult = findNextEdge();
			currentEdge = nextEdgeResult.getKey();
			boolean endBranchMic = nextEdgeResult.getValue();

			if (currentEdge == null) {
				return false;
			}

			if (endBranchMic) {
				return true;
			}

			Entry<Double, Double> nextUR = findNextU();
			if (newBranch) {
				newBranch = false;
				outputNextMic(nextUR.getKey(), nextUR.getValue(), true);
			} else {
				outputNextMic(nextUR.getKey(), nextUR.getValue(), false);
			}
			return true;
		}
	}

	/**
	 * Pop an unvisited edge from the stack, or return {@code null} if the stack
	 * is empty (end-of-operation).
	 */
	private Edge findNextBranch() {
		if (unvisited.isEmpty()) {
			return null;
		}
		BranchPoint bp = unvisited.pop();
		previousBranchCenter = currentCenter;
		previousBranchRadius = currentRadius;
		currentCenter = bp.currentCenter;
		currentRadius = bp.currentRadius;
		newBranch = true;
		if (!edgeDone.getOrDefault(bp.nextEdge, false)) {
			return bp.nextEdge;
		} else {
			return findNextBranch();
		}
	}

	/** Find valid out-edges of the current edge's target vertex. */
	private List<Edge> findOutEdges() {
		Vertex trg = currentEdge.target;
		List<Edge> outEdges = new ArrayList<>();
		for (Edge e : trg.outEdges) {
			if (e != currentEdge.twin
					&& e.valid
					&& !edgeDone.getOrDefault(e, false)
					&& e.type != EdgeType.NULLEDGE
					&& e.type != EdgeType.OUTEDGE) {
				outEdges.add(e);
			}
		}
		return outEdges;
	}

	/**
	 * Find the next edge to advance to. Returns a pair of (edge, endBranchMic).
	 * If edge is {@code null}, the operation is over. If endBranchMic is
	 * {@code true}, a final MIC at the end of a branch is needed.
	 */
	private Entry<Edge, Boolean> findNextEdge() {
		List<Edge> outEdges = findOutEdges();

		if (outEdges.isEmpty()) {
			// end of branch — possibly output a final MIC
			if (currentRadius > currentEdge.target.dist()) {
				currentRadius = currentEdge.target.dist();
				currentCenter = currentEdge.target.position;
				return new SimpleEntry<>(currentEdge, true);
			}
			Edge e = findNextBranch();
			if (e == null) {
				return new SimpleEntry<>(null, false);
			}
			if (hasNextRadius(e)) {
				currentU = 0;
				return new SimpleEntry<>(e, false);
			} else {
				markDone(e);
				currentEdge = e;
				return findNextEdge();
			}
		} else if (outEdges.size() == 1) {
			if (hasNextRadius(outEdges.get(0))) {
				currentU = 0;
				return new SimpleEntry<>(outEdges.get(0), false);
			} else {
				markDone(outEdges.get(0));
				currentEdge = outEdges.get(0);
				return findNextEdge();
			}
		} else if (outEdges.size() == 2) {
			unvisited.push(new BranchPoint(currentCenter, currentRadius, outEdges.get(1)));
			if (hasNextRadius(outEdges.get(0))) {
				currentU = 0;
				return new SimpleEntry<>(outEdges.get(0), false);
			} else {
				markDone(outEdges.get(0));
				currentEdge = outEdges.get(0);
				return findNextEdge();
			}
		} else {
			// More than 2 out-edges — push extras, take the first
			for (int i = 1; i < outEdges.size(); i++) {
				unvisited.push(new BranchPoint(currentCenter, currentRadius, outEdges.get(i)));
			}
			if (hasNextRadius(outEdges.get(0))) {
				currentU = 0;
				return new SimpleEntry<>(outEdges.get(0), false);
			} else {
				markDone(outEdges.get(0));
				currentEdge = outEdges.get(0);
				return findNextEdge();
			}
		}
	}

	/** Mark an edge and its twin as done. */
	private void markDone(Edge e) {
		edgeDone.put(e, true);
		if (e.twin != null) {
			edgeDone.put(e.twin, true);
		}
	}

	/**
	 * Does edge {@code e} have a next MIC at its target that exceeds the
	 * cut-width threshold?
	 */
	private boolean hasNextRadius(Edge e) {
		Entry<Point, Double> sample = edgePoint(e, 1.0);
		double w = cutWidth(currentCenter, currentRadius, sample.getKey(), sample.getValue());
		if (w <= 0) {
			return false;
		}
		return w > maxWidth;
	}

	/**
	 * Find the next u-value on the current edge that produces the desired
	 * cut-width. Uses bracketed root-finding (Brent's method).
	 *
	 * @return pair of (nextU, nextRadius)
	 */
	private Entry<Double, Double> findNextU() {
		CutWidthError errFunc = new CutWidthError(currentEdge, maxWidth, currentCenter, currentRadius);

		BracketingNthOrderBrentSolver solver = new BracketingNthOrderBrentSolver(1e-9, 5);
		double result = solver.solve(500, errFunc, currentU, 1.0, AllowedSolution.ANY_SIDE);

		Entry<Point, Double> sample = edgePoint(currentEdge, result);
		return new SimpleEntry<>(result, sample.getValue());
	}

	/**
	 * Output the next MIC and advance the current state.
	 */
	private void outputNextMic(double nextU, double nextRadius, boolean branch) {
		MIC mic = new MIC();
		Point c1 = currentCenter;
		double r1 = currentRadius;
		Entry<Point, Double> sample = edgePoint(currentEdge, nextU);
		Point c2 = sample.getKey();
		double r2 = sample.getValue();

		mic.c1 = c1;
		mic.r1 = r1;
		mic.c2 = c2;
		mic.r2 = r2;

		if (!c1.equals(c2)) {
			List<Point> tangents = bitangentPoints(c1, r1, c2, r2);
			mic.t1 = tangents.get(0);
			mic.t2 = tangents.get(1);
			mic.t3 = tangents.get(2);
			mic.t4 = tangents.get(3);
		}
		mic.newBranch = branch;
		mic.cPrev = previousBranchCenter;
		mic.rPrev = previousBranchRadius;
		micList.add(mic);

		currentRadius = r2;
		currentCenter = c2;
		currentU = nextU;
		
		Point[] fp1 = edgeFootPoints(currentEdge, currentU);
		if (fp1 != null && fp1.length == 2) {
			mic.tp1a = fp1[0];
			mic.tp1b = fp1[1];
		}

		Point[] fp2 = edgeFootPoints(currentEdge, nextU);
		if (fp2 != null && fp2.length == 2) {
			mic.tp2a = fp2[0];
			mic.tp2b = fp2[1];
		}
	}

	/**
	 * Compute four bi-tangent points between two circles (c1,r1) and (c2,r2).
	 */
	static List<Point> bitangentPoints(Point c1, double r1, Point c2, double r2) {
		List<Point> out = new ArrayList<>(4);
		Point bd1, bd2;

		if (r1 == r2) {
			Point c1c2 = c2.sub(c1);
			c1c2.normalize();
			bd1 = c1c2.xyPerp().mult(-1);
			bd2 = c1c2.xyPerp();
		} else {
			double c1c2Dist = c1.sub(c2).norm();
			double dr = Math.abs(r1 - r2);
			double bitangLength = Math.sqrt(c1c2Dist * c1c2Dist + dr * dr);
			double area = dr * bitangLength;
			double height = area / c1c2Dist;
			double bitangC1c2 = Math.sqrt(bitangLength * bitangLength - height * height);
			Point cdir;
			if (r1 > r2) {
				cdir = c1.sub(c2);
			} else {
				cdir = c2.sub(c1);
			}
			cdir.normalize();

			Point bit1 = cdir.mult(bitangC1c2).add(cdir.xyPerp().mult(height));
			Point bit2 = cdir.mult(bitangC1c2).sub(cdir.xyPerp().mult(height));

			if (r1 > r2) {
				bd1 = bit1.xyPerp();
				bd2 = bit2.xyPerp().mult(-1);
			} else {
				bd1 = bit2.xyPerp().mult(-1);
				bd2 = bit1.xyPerp();
			}
			bd1.normalize();
			bd2.normalize();
		}

		out.add(c1.add(bd1.mult(r1)));
		out.add(c1.add(bd2.mult(r1)));
		out.add(c2.add(bd1.mult(r2)));
		out.add(c2.add(bd2.mult(r2)));
		return out;
	}

	/**
	 * Return cut-width from (c1,r1) to (c2,r2).
	 * <p>
	 * {@code (c1, r1)} is the previously machined MIC and {@code (c2, r2)} is
	 * the new MIC. The maximum cut-width when cutting c2 is:
	 * <pre>
	 *   w_max = |c2 - c1| + r2 - r1
	 * </pre>
	 */
	static double cutWidth(Point c1, double r1, Point c2, double r2) {
		return c2.sub(c1).norm() + r2 - r1;
	}

	/**
	 * Find a point on the given edge at normalized parameter {@code u ∈ [0,1]}.
	 * This delegates to the existing {@link Edge#micSample(double)} method.
	 *
	 * @param e edge to sample
	 * @param u normalized parameter from 0.0 (source) to 1.0 (target)
	 * @return pair of (point, clearanceRadius)
	 */
	private static Entry<Point, Double> edgePoint(Edge e, double u) {
		return e.micSample(u);
	}
	
	private Point[] edgeFootPoints(Edge e, double u) {
		return e.micTouchPoints(u);
	}

	/**
	 * Error functor for {@link #findNextU()}: returns the difference between
	 * the cut-width at parameter {@code x} and the desired maximum cut-width.
	 */
	private class CutWidthError implements UnivariateFunction {
		private final Edge e;
		private final double wMax;
		private final Point c1;
		private final double r1;

		CutWidthError(Edge e, double wMax, Point c1, double r1) {
			this.e = e;
			this.wMax = wMax;
			this.c1 = c1;
			this.r1 = r1;
		}

		@Override
		public double value(double x) {
			Entry<Point, Double> sample = edgePoint(e, x);
			Point c2 = sample.getKey();
			double r2 = sample.getValue();
			double w = c2.sub(c1).norm() + r2 - r1;
			return w - wMax;
		}
	}

	/**
	 * Branch-data when backtracking to machine an un-machined branch.
	 */
	private static class BranchPoint {
		final Point currentCenter;
		final double currentRadius;
		final Edge nextEdge;

		BranchPoint(Point p, double r, Edge e) {
			this.currentCenter = p;
			this.currentRadius = r;
			this.nextEdge = e;
		}
	}
}
