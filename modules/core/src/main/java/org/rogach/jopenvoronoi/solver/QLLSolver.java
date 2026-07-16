package org.rogach.jopenvoronoi.solver;

import static org.rogach.jopenvoronoi.util.Numeric.chop;
import static org.rogach.jopenvoronoi.util.Numeric.quadraticRoots;

import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Quadratic-linear-linear solver.
 */
public class QLLSolver extends Solver {

	// reusable per-solve scratch; solver instances are used single-threaded
	// (see VertexPositioner.solutionsBuffer for the same reasoning)
	private final Eq[] eqScratch = { new Eq(), new Eq(), new Eq() };
	private final Eq[] quads = new Eq[3];
	private final Eq[] lins = new Eq[3];
	// per-qllSolver-call scratch; every element is written before it is read
	private final double[][] aargs = new double[3][2];
	private final double[][] isolns = new double[2][3];
	private final double[][] tsolns = new double[2][3];
	private final double[] roots = new double[2];

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		// equation-parameters, in quad-precision
		var nQuads = 0;
		var nLins = 0;
		for (var i = 0; i < 3; i++) {
			var site = (i == 0) ? s1 : (i == 1) ? s2 : s3;
			var kval = (i == 0) ? k1 : (i == 1) ? k2 : k3;
			var eqn = eqScratch[i];
			site.eqp(kval, eqn);
			if (site.isLine()) {
				lins[nLins++] = eqn;
			} else {
				quads[nQuads++] = eqn;
			}
		}
		assert (nQuads > 0) : " !quads.isEmpty() ";

		if (nLins == 1 || nLins == 0) {
			assert (nQuads == 3 || nQuads == 2) : " quads.size() == 3 || quads.size() == 2 ";
			for (var i = 1; i < nQuads; i++) {
				quads[i].subEq(quads[0]);
				lins[nLins++] = quads[i];
			}
		}
		assert (nLins == 2) : " lins.size() == 2";

		// TODO: pick the solution appraoch with the best numerical stability.
		// call all three permutations
		// index shuffling determines if we solve:
		// x and y in terms of t
		// y and t in terms of x
		// t and x in terms of y
		qllSolver(lins, 0, 1, 2, quads[0], k3, slns);
		qllSolver(lins, 2, 0, 1, quads[0], k3, slns);
		qllSolver(lins, 1, 2, 0, quads[0], k3, slns);

		return slns.size();
	}

	/**
	 * Solve the QLL system for a fixed quadratic equation and two linear equations.
	 *
	 * @param lins two linear equations
	 * @param xi x-index used when shuffling coordinates
	 * @param yi y-index used when shuffling coordinates
	 * @param ti t-index used when shuffling coordinates
	 * @param quad parameters of the quadratic site (point or arc)
	 * @param k3 offset direction for the quadratic site
	 * @param solns output solution triplets {@code (x, y, t)} or
	 *              {@code (u, v, t)}
	 * @return number of solutions found
	 */
	private int qllSolver(Eq[] lins, int xi, int yi, int ti, Eq quad, double k3, List<Solution> solns) {
		var ai = lins[0].get(xi); // first linear
		var bi = lins[0].get(yi);
		var ki = lins[0].get(ti);
		var ci = lins[0].c;

		var aj = lins[1].get(xi); // second linear
		var bj = lins[1].get(yi);
		var kj = lins[1].get(ti);
		var cj = lins[1].c;

		var d = chop(ai * bj - aj * bi); // chop! (determinant for 2 linear eqns (?))
		if (d == 0) {
			return -1;
		}
		// these are the w-equations for qll_solve()
		// (2) u = a1 w + b1
		// (3) v = a2 w + b2
		var a0 = (bi * kj - bj * ki) / d;
		var a1 = -(ai * kj - aj * ki) / d;
		var b0 = (bi * cj - bj * ci) / d;
		var b1 = -(ai * cj - aj * ci) / d;
		// based on the 'last' quadratic of (s1,s2,s3)
		aargs[0][0] = 1.0;
		aargs[0][1] = quad.a;
		aargs[1][0] = 1.0;
		aargs[1][1] = quad.b;
		aargs[2][0] = -1.0;
		aargs[2][1] = quad.k;

		// this solves for w, and returns either 0, 1, or 2 triplets of (u,v,t) in
		// isolns
		// NOTE: indexes of aargs shuffled depending on (xi,yi,ti) !
		var scount = qllSolve(aargs[xi][0], aargs[xi][1], aargs[yi][0], aargs[yi][1], aargs[ti][0], aargs[ti][1],
				quad.c, // xk*xk + yk*yk - rk*rk,
				a0, b0, a1, b1, isolns);
		for (var i = 0; i < scount; i++) {
			tsolns[i][xi] = isolns[i][0]; // u x
			tsolns[i][yi] = isolns[i][1]; // v y
			tsolns[i][ti] = isolns[i][2]; // t t chop!
			solns.add(new Solution(new Point(tsolns[i][0], tsolns[i][1]), tsolns[i][2], k3));
		}
		// std::cout << " k3="<<kk3<<" qqq_solve found " << scount << " roots\n";
		return scount;
	}

	// Solve a system of one quadratic equation, and two linear equations.
	/**
	 * Solve a system of one quadratic equation and two linear equations:
	 * <pre>{@code
	 * (1) a0 u^2 + b0 u + c0 v^2 + d0 v + e0 w^2 + f0 w + g0 = 0
	 * (2) u = a1 w + b1
	 * (3) v = a2 w + b2
	 * }</pre>
	 * Solve equation (1) for {@code w} and then substitute into (2) and (3) to
	 * find {@code (u, v, t)}.
	 */
	private int qllSolve(double a0, double b0, double c0, double d0, double e0, double f0, double g0, double a1,
			double b1, double a2, double b2, double soln[][]) {
		// std::cout << "qll_solver()\n";
		// TODO: optimize using abs(a0) == abs(c0) == abs(d0) == 1
		var a = chop((a0 * (a1 * a1) + c0 * (a2 * a2) + e0));
		var b = chop((2 * a0 * a1 * b1 + 2 * a2 * b2 * c0 + a1 * b0 + a2 * d0 + f0));
		var c = a0 * (b1 * b1) + c0 * (b2 * b2) + b0 * b1 + b2 * d0 + g0;
		var nRoots = quadraticRoots(a, b, c, roots); // solves a*w^2 + b*w + c = 0
		for (var i = 0; i < nRoots; i++) {
			double w = roots[i];
			soln[i][0] = a1 * w + b1; // u
			soln[i][1] = a2 * w + b2; // v
			soln[i][2] = w; // t
		}
		return nRoots;
	}

}
