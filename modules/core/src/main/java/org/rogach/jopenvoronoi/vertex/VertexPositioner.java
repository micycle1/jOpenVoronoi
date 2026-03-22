package org.rogach.jopenvoronoi.vertex;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.solver.ALTSEPSolver;
import org.rogach.jopenvoronoi.solver.LLLPARASolver;
import org.rogach.jopenvoronoi.solver.LLLSolver;
import org.rogach.jopenvoronoi.solver.PPPSolver;
import org.rogach.jopenvoronoi.solver.QLLSolver;
import org.rogach.jopenvoronoi.solver.SEPSolver;
import org.rogach.jopenvoronoi.solver.Solver;
import org.rogach.jopenvoronoi.util.Numeric;

/**
 * Calculates the {@code (x, y)} position of a Voronoi vertex in the
 * {@link org.rogach.jopenvoronoi.VoronoiDiagram Voronoi diagram}.
 */
public class VertexPositioner {

	// solvers, to which we dispatch, depending on the input sites
	/** point-point-point solver */
	Solver pppSolver;
	/** line-line-line solver */
	Solver lllSolver;
	/** parallel line-line-line solver */
	Solver lllParaSolver;
	/** quadratic-linear-linear solver */
	Solver qllSolver;
	/** separator solver */
	Solver sepSolver;
	/** alternative separator solver */
	Solver altSepSolver;

	// DATA
	/** reference to the VD graph. */
	HalfEdgeDiagram g;
	/** minimum offset-distance */
	double tMin;
	/** maximum offset-distance */
	double tMax;
	/** the edge on which we position a new vertex */
	Edge edge;

	// create positioner, set graph.
	public VertexPositioner(HalfEdgeDiagram gi) {
		this.g = gi;
		pppSolver = new PPPSolver();
		lllSolver = new LLLSolver();
		qllSolver = new QLLSolver();
		sepSolver = new SEPSolver();
		altSepSolver = new ALTSEPSolver();
		lllParaSolver = new LLLPARASolver();
	}

	/**
	 * Position a new vertex on the given half-edge {@code e} when inserting the new
	 * Site {@code s3}.
	 * <p>
	 * The new vertex is equidistant to the two sites that defined the edge and to
	 * the new site. The edge {@code e} holds information about which face it
	 * belongs to, and each face holds information about which site created it.
	 * <p>
	 * The three sites defining the position of the vertex are:
	 * <ul>
	 * <li>the site to the left of half-edge {@code e}</li>
	 * <li>the site to the right of half-edge {@code e}</li>
	 * <li>the given new Site {@code s3}</li>
	 * </ul>
	 */
	public Solution position(Edge e, Site s3) {
		edge = e;
		var face = e.face;
		var twin = e.twin;
		var twin_face = twin.face;

		var src = e.source;
		var trg = e.target;
		var t_src = src.dist();
		var t_trg = trg.dist();
		tMin = Math.min(t_src, t_trg); // the solution we seek must have t_min<t<t_max
		tMax = Math.max(t_src, t_trg);

		var s1 = face.getSite();
		var s2 = twin_face.getSite();

		var sl = position(s1, e.k, s2, twin.k, s3);

		assert (solution_on_edge(sl)) : " solution_on_edge(sl) ";
		// assert( check_far_circle(sl) ) : " check_far_circle(sl) ";
		assert (check_dist(edge, sl, s3)) : " check_dist(edge, sl, s3) ";

		return sl;
	}

	// position new vertex
	// find vertex that is equidistant from s1, s2, s3
	// should lie on the k1 side of s1, k2 side of s2
	// we try both k3=-1 and k3=+1 for s3
	Solution position(Site s1, double k1, Site s2, double k2, Site s3) {
		assert ((k1 == 1) || (k1 == -1)) : " (k1==1) || (k1 == -1) ";
		assert ((k2 == 1) || (k2 == -1)) : " (k2==1) || (k2 == -1) ";

		if (s3.isLine()) {
			// special handling for the case when site and edge endpoints share a common
			// point -
			// simply select one of the edge endpoints as solution
			var e_src = edge.source.position;
			if (((s1.isPoint() && s1.position().equals(e_src))
					|| (s1.isLine() && (s1.start().equals(e_src) || s1.end().equals(e_src))))
					&& ((s2.isPoint() && s2.position().equals(e_src))
							|| (s2.isLine() && (s2.start().equals(e_src) || s2.end().equals(e_src))))
					&& ((s3.isPoint() && s3.position().equals(e_src))
							|| (s3.isLine() && (s3.start().equals(e_src) || s3.end().equals(e_src))))) {
				var src_se = s3.start();
				var trg_se = s3.end();
				double k;
				if (edge.target.position.isRight(src_se, trg_se)) {
					k = (s3.k() == 1) ? -1 : 1;
				} else {
					k = (s3.k() == 1) ? 1 : -1;
				}
				return new Solution(edge.source.position, edge.source.dist(), k);
			}
			var e_trg = edge.target.position;
			if (((s1.isPoint() && s1.position().equals(e_trg))
					|| (s1.isLine() && (s1.start().equals(e_trg) || s1.end().equals(e_trg))))
					&& ((s2.isPoint() && s2.position().equals(e_trg))
							|| (s2.isLine() && (s2.start().equals(e_trg) || s2.end().equals(e_trg))))
					&& ((s3.isPoint() && s3.position().equals(e_trg))
							|| (s3.isLine() && (s3.start().equals(e_trg) || s3.end().equals(e_trg))))) {
				var src_se = s3.start();
				var trg_se = s3.end();
				double k;
				if (edge.source.position.isRight(src_se, trg_se)) {
					k = (s3.k() == 1) ? -1 : 1;
				} else {
					k = (s3.k() == 1) ? 1 : -1;
				}
				return new Solution(edge.target.position, edge.target.dist(), k);
			}
		}

		List<Solution> solutions = new ArrayList<>();

		if (s3.isLine() && ((s1.isPoint() && s2.isLine()
				&& (s3.start().equals(s1.position()) || s3.end().equals(s1.position())))
				|| (s2.isPoint() && s1.isLine()
						&& (s3.start().equals(s2.position()) || s3.end().equals(s2.position()))))) {
			var ptsite = s1.isPoint() ? s1 : s2;
			var ed = edge;
			if (ed.face.getSite() != ptsite) {
				ed = edge.twin;
			}
			assert (ed.source.status == VertexStatus.IN || ed.target.status == VertexStatus.IN)
					: "edge to be split has no IN vertex";

			double k;
			if (ed.source.status == VertexStatus.IN) {
				k = -1;
			} else {
				k = +1;
			}
			if (s3.start().equals(ptsite.position())) {
				k = -k;
			}
			solverDispatch(s1, k1, s2, k2, s3, k, solutions);
		} else {
			solverDispatch(s1, k1, s2, k2, s3, +1, solutions); // a single k3=+1 call for s3->isPoint()

			if (!s3.isPoint()) {
				solverDispatch(s1, k1, s2, k2, s3, -1, solutions); // for lineSite or ArcSite we try k3=-1 also
			}
		}

		if (solutions.size() == 1 && (tMin <= solutions.get(0).t) && (tMax >= solutions.get(0).t)
				&& (s3.inRegion(solutions.get(0).p))) {
			return solutions.get(0);
		}

		// choose only in_region() solutions
		List<Solution> acceptable_solutions = new ArrayList<>();
		for (Solution s : solutions) {
			if (s3.inRegion(s.p) && s.t >= tMin && s.t <= tMax) {
				acceptable_solutions.add(s);
			}
		}

		if (acceptable_solutions.size() == 1) { // if only one solution is found, return that.
			return acceptable_solutions.get(0);
		} else if (acceptable_solutions.size() > 1) {
			// two or more points remain so we must further filter here!
			// filter further using edge_error
			var min_error = 100D;
			var min_solution = new Solution(new Point(0, 0), 0, 0);
			for (Solution s : acceptable_solutions) {
				var err = edgeError(s);
				if (err < min_error) {
					min_solution = s;
					min_error = err;
				}
			}
			return min_solution;
		}

		if (solutions.isEmpty()) {
			return desperateSolution(s3);
		} else {
			// choose solution that is best by dist_error
			var leastBad = solutions.get(0);
			var leastErr = Double.MAX_VALUE;
			for (Solution s : solutions) {
				// punish wrong solutions
				var derr = distError(edge, s, s3);
				// punish solutions outside t range
				var terr = Math.max(0, Math.max((s.t - tMax), (tMin - s.t)));
				if (edge.type == EdgeType.PARA_LINELINE) {
					var s_p = s.p.sub(edge.source.position);
					var s_e = edge.target.position.sub(edge.source.position);
					var dist = s_p.dot(s_e) / s_e.dot(s_e);
					terr = Math.max(0, Math.max(dist - 1, -dist));
				}
				var err = derr + terr;
				if (err < leastErr) {
					leastBad = s;
					leastErr = err;
				}
			}

			if (edge.type == EdgeType.PARA_LINELINE) {
				return leastBad;
			}

			// determine clamp direction
			var t = Math.max(tMin, Math.min(tMax, leastBad.t));
			var p_sln = edge.point(t);

			// find out on which side the solution lies
			var desp_k3 = 0D;
			if (s3.isPoint()) {
				desp_k3 = 1;
			} else if (s3.isLine()) {
				var src_se = s3.start();
				var trg_se = s3.end();
				if (p_sln.isRight(src_se, trg_se)) {
					desp_k3 = (s3.k() == 1) ? -1 : 1;
				} else {
					desp_k3 = (s3.k() == 1) ? 1 : -1;
				}
			}
			return new Solution(p_sln, t, desp_k3);
		}
	}

	// search numerically for a desperate solution along the solution-edge
	Solution desperateSolution(Site s3) {
		var err_functor = new VertexError(g, edge, s3);

		double t_sln;
		// Guard against zero-length or negative search intervals
		if (Math.abs(tMax - tMin) < Numeric.STRICT_ZERO_EPSILON) {
			t_sln = tMin; // Interval is effectively zero, skip the optimizer
		} else {
			var optimizer = new BrentOptimizer(Numeric.DEFAULT_CHOP_EPSILON, Numeric.STRICT_ZERO_EPSILON);
			t_sln = optimizer.optimize(
					new MaxEval(1000), 
					new UnivariateObjectiveFunction(err_functor),
					GoalType.MINIMIZE, 
					new SearchInterval(tMin, tMax)
			).getPoint();
		}
		var p_sln = err_functor.edgePoint(t_sln); // g[edge].point(t_sln);
		var desp_k3 = 0D;
		if (s3.isPoint()) {
			desp_k3 = 1;
		} else if (s3.isLine()) {
			// find out on which side the desperate solution lies
			var src_se = s3.start();
			var trg_se = s3.end();
			if (p_sln.isRight(src_se, trg_se)) {
				desp_k3 = (s3.k() == 1) ? -1 : 1;
			} else {
				desp_k3 = (s3.k() == 1) ? 1 : -1;
			}
		}
		return new Solution(p_sln, t_sln, desp_k3);
	}

	// dispatch to the correct solver based on the sites
	int solverDispatch(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> solns) {

		if (edge.type == EdgeType.SEPARATOR) {
			// this is a SEPARATOR edge with two LineSites adjacent.
			// find the PointSite that defines the SEPARATOR, so that one LineSite and one
			// PointSite
			// can be submitted to the Solver.
			if (s1.isLine() && s2.isLine()) {
				// the parallell lineseg case v0 --s1 --> pt -- s2 --> v1
				// find t
				if (edge.hasNullFace) {
					s2 = edge.nullFace.getSite();
					assert (s2.isPoint()) : " s2.isPoint() ";
					k2 = +1;
				} else if (edge.twin.hasNullFace) {
					s2 = edge.twin.nullFace.getSite();
					assert (s2.isPoint()) : " s2.isPoint() ";
					k2 = +1;
				}
			} else if (s1.isPoint() && s2.isLine()) {
				// a normal SEPARATOR edge, defined by a PointSite and a LineSite
				// swap sites, so SEPSolver can assume s1=line s2=point
				var tmp = s1;
				var k_tmp = k1;
				s1 = s2;
				s2 = tmp;
				k1 = k2;
				k2 = k_tmp;
				assert (s1.isLine()) : " s1.isLine() ";
				assert (s2.isPoint()) : " s2.isPoint() ";
			}
			assert (s1.isLine() && s2.isPoint()) : " s1.isLine() && s2.isPoint() ";
			return sepSolver.solve(s1, k1, s2, k2, s3, k3, solns);
		} else if (edge.type == EdgeType.PARA_LINELINE && s3.isLine()) { // an edge betwee parallel LineSites
			// std::cout << " para lineline! \n";
			return lllParaSolver.solve(s1, k1, s2, k2, s3, k3, solns);
		} else if (s1.isLine() && s2.isLine() && s3.isLine()) {
			return lllSolver.solve(s1, k1, s2, k2, s3, k3, solns); // all lines.
		} else if (s1.isPoint() && s2.isPoint() && s3.isPoint()) {
			return pppSolver.solve(s1, 1, s2, 1, s3, 1, solns); // all points, no need to specify k1,k2,k3, they are
																	// all +1
		} else if ((s3.isLine() && s1.isPoint()) || (s1.isLine() && s3.isPoint()) || (s3.isLine() && s2.isPoint())
				|| (s2.isLine() && s3.isPoint()) // bad coverage for this line?
		) {
			// if s1/s2 form a SEPARATOR-edge, this is dispatched automatically to
			// sep-solver
			// here we detect for a separator case between
			// s1/s3
			// s2/s3
			if (s3.isLine() && s1.isPoint()) {
				if (detectSepCase(s3, s1)) {
					altSepSolver.setType(0);
					return altSepSolver.solve(s1, k1, s2, k2, s3, k3, solns);
				}
			}
			if (s3.isLine() && s2.isPoint()) {
				if (detectSepCase(s3, s2)) {
					altSepSolver.setType(1);
					return altSepSolver.solve(s1, k1, s2, k2, s3, k3, solns);
				}
			}
		}

		// if we didn't dispatch to a solver above, we try the general solver
		return qllSolver.solve(s1, k1, s2, k2, s3, k3, solns); // general case solver

	}

	// detect separator-case, so we can dispatch to the correct Solver
	boolean detectSepCase(Site lsite, Site psite) {
		var le = lsite.edge();
		var src = le.source;
		var trg = le.target;
		// now from segment end-points get the null-vertex
		Edge src_out = null;
		for (Edge e : src.outEdges) {
			if (e.type == EdgeType.NULLEDGE) {
				src_out = e;
			}
		}
		Edge trg_out = null;
		for (Edge e : trg.outEdges) {
			if (e.type == EdgeType.NULLEDGE) {
				trg_out = e;
			}
		}

		var src_null_face = src_out.face;
		if (src_null_face.isNullFace() == false) {
			// take twin face instead
			var src_out_twin = src_out.twin;
			src_null_face = src_out_twin.face;
		}

		var trg_null_face = trg_out.face;
		if (trg_null_face.isNullFace() == false) {
			var trg_out_twin = trg_out.twin;
			trg_null_face = trg_out_twin.face;
		}
		assert (src_null_face.isNullFace() && trg_null_face.isNullFace())
				: " src_null_face.isNullFace() && trg_null_face.isNullFace() ";

		// do we want src_out face??
		// OR src_out_twin face??
		// we want the null-face !

		var src_site = src_null_face.getSite();
		var trg_site = trg_null_face.getSite();
		if (src_site == null || trg_site == null) {
			throw new RuntimeException();
		}
		if (!src_site.isPoint() || !trg_site.isPoint()) {
			throw new RuntimeException();
		}
		var src_vertex = src_site.vertex();
		var trg_vertex = trg_site.vertex();
		if (src_vertex == psite.vertex()) {
			return true;
		}
		if (trg_vertex == psite.vertex()) {
			return true;
		}
		return false;
	}

	// error from solution to corresponding point on the edge
	double edgeError(Solution sl) {
		Point p;
		if (edge.type == EdgeType.PARA_LINELINE) {
			p = projectionPoint(sl);
		} else {
			p = edge.point(sl.t);
		}
		return p.sub(sl.p).norm();
	}

	// when the edge is not parametrized by t-value as normal edges
	// so we need a projection of sl onto the edge instead
	Point projectionPoint(Solution sl) {
		assert (edge.type == EdgeType.PARA_LINELINE) : " edge.type == EdgeType.PARA_LINELINE ";
		// edge given by
		// p = p0 + t * (p1-p0) with t in [0,1]
		var p0 = new Point(edge.source.position);
		var p1 = new Point(edge.target.position);
		var v = p1.sub(p0);

		var t = sl.p.sub(p0).dot(v) / v.dot(v);
		// clamp to [0,1]
		if (t > 1) {
			t = 1;
		} else if (t < 0) {
			t = 0;
		}
		return p0.add(v.mult(t));
	}

	// check that the new solution lies on the edge
	boolean solution_on_edge(Solution s) {
		var err = edgeError(s);
		var limit = Numeric.SOLUTION_EDGE_EPSILON;
		return (err < limit);
	}

	// new vertices should lie within the far_radius
	boolean checkFarCircle(Solution s) {
		if (!(s.p.norm() < 18 * 1)) {
			return false;
		}
		return true;
	}

	// distance sanity check
	// all vertices should be of degree three, i.e. three adjacent faces/sites
	// distance to the three adjacent sites should be equal
	boolean check_dist(Edge e, Solution sl, Site s3) {
		var face = e.face;
		var tw_edge = e.twin;
		var twin_face = tw_edge.face;

		var s1 = face.getSite();
		var s2 = twin_face.getSite();

		var d1 = sl.p.sub(s1.apexPoint(sl.p)).norm();
		var d2 = sl.p.sub(s2.apexPoint(sl.p)).norm();
		var d3 = sl.p.sub(s3.apexPoint(sl.p)).norm();

		if (!equal(d1, d2) || !equal(d1, d3) || !equal(d2, d3) || !equal(sl.t, d1) || !equal(sl.t, d2)
				|| !equal(sl.t, d3)) {
			return false;
		}
		return true;
	}

	// distance-error
	// new vertices should be equidistant to the three adjacent sites that define
	// the vertex
	// we here calculate the distances d1, d2, d3 from the Solution to the three
	// sites s1, s2, s3
	// and return the max deviation from the solution t-value.
	// this works as a sanity check for the solver.
	// a high error value here is also an indication of numerical instability in the
	// solver
	public double distError(Edge e, Solution sl, Site s3) {
		var face = e.face;
		var tw_edge = e.twin;
		var twin_face = tw_edge.face;

		var s1 = face.getSite();
		var s2 = twin_face.getSite();

		var d1 = sl.p.sub(s1.apexPoint(sl.p)).norm();
		var d2 = sl.p.sub(s2.apexPoint(sl.p)).norm();
		var d3 = sl.p.sub(s3.apexPoint(sl.p)).norm();

		return Math.max(Math.max(Math.abs(sl.t - d1), Math.abs(sl.t - d2)), Math.abs(sl.t - d3));

	}

	// are \a d1 and \a d2 roughly equal?
	boolean equal(double d1, double d2) {
		return Numeric.areClose(d1, d2, Numeric.DISTANCE_EPSILON, Numeric.DOUBLE_COMPARISON_EPSILON);
	}

}
