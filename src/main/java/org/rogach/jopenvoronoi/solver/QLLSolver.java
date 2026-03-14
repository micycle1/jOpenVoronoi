package org.rogach.jopenvoronoi.solver;

import java.util.ArrayList;
import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Quadratic-linear-linear solver.
 * <p>
 * Handles one quadratic equation and two linear equations. If there are
 * multiple quadratic sites, we subtract one quadratic equation from the others
 * to obtain the needed linear equations.
 */
public class QLLSolver extends Solver {

	private static final int[][] PERMUTATIONS = { { 0, 1, 2 }, // solve x,y in terms of t
			{ 2, 0, 1 }, // solve t,x in terms of y
			{ 1, 2, 0 } // solve y,t in terms of x
	};

	private static final double DET_EPS = 1e-12;
	private static final double ROOT_EPS = 1e-12;
	private static final double RES_EPS = 1e-9;
	private static final double DUP_EPS = 1e-9;

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		var sites = new Site[] { s1, s2, s3 };
		var kvals = new double[] { k1, k2, k3 };

		List<Eq> original = new ArrayList<>(3);
		List<Eq> quads = new ArrayList<>(3);
		List<Eq> lins = new ArrayList<>(2);

		for (int i = 0; i < 3; i++) {
			var eq = sites[i].eqp(kvals[i]);
			original.add(new Eq(eq));
			if (sites[i].is_linear()) {
				lins.add(new Eq(eq));
			} else {
				quads.add(new Eq(eq));
			}
		}

		if (quads.isEmpty()) {
			return 0;
		}

		// Convert extra quadratic equations into linear equations by subtraction.
		if (lins.size() <= 1) {
			if (quads.size() < 2) {
				return 0;
			}
			var base = quads.get(0);
			for (int i = 1; i < quads.size(); i++) {
				var lin = new Eq(quads.get(i));
				lin.subEq(base);
				lins.add(lin);
			}
		}

		if (lins.size() != 2) {
			return 0;
		}

		var quad = quads.get(0);

		// Try best-conditioned permutation first, then the others if needed.
		var order = permutationOrder(lins);
		for (int idx : order) {
			var p = PERMUTATIONS[idx];
			int added = solvePermutation(lins, p[0], p[1], p[2], quad, original, k3, slns);
			if (added > 0) {
				return added;
			}
		}

		return 0;
	}

	private int solvePermutation(List<Eq> lins, int xi, int yi, int ti, Eq quad, List<Eq> original, double k3, List<Solution> solns) {

		var li = lins.get(0);
		var lj = lins.get(1);

		double ai = li.get(xi);
		double bi = li.get(yi);
		double ki = li.get(ti);
		double ci = li.c;

		double aj = lj.get(xi);
		double bj = lj.get(yi);
		double kj = lj.get(ti);
		double cj = lj.c;

		double d = ai * bj - aj * bi;
		double dScale = Math.abs(ai * bj) + Math.abs(aj * bi) + 1.0;
		if (Math.abs(d) <= DET_EPS * dScale) {
			return 0;
		}

		// u = uSlope*w + uIntercept
		// v = vSlope*w + vIntercept
		double uSlope = (bi * kj - bj * ki) / d;
		double vSlope = -(ai * kj - aj * ki) / d;
		double uIntercept = (bi * cj - bj * ci) / d;
		double vIntercept = -(ai * cj - aj * ci) / d;

		double[] squareCoeff = { 1.0, 1.0, -1.0 };
		double[] linearCoeff = { quad.a, quad.b, quad.k };

		double[] roots = new double[2];
		int rootCount = solveReducedQuadratic(squareCoeff[xi], linearCoeff[xi], squareCoeff[yi], linearCoeff[yi], squareCoeff[ti], linearCoeff[ti], quad.c,
				uSlope, uIntercept, vSlope, vIntercept, roots);

		int added = 0;
		for (int i = 0; i < rootCount; i++) {
			double w = roots[i];
			double u = uSlope * w + uIntercept;
			double v = vSlope * w + vIntercept;

			double[] xyz = new double[3];
			xyz[xi] = u;
			xyz[yi] = v;
			xyz[ti] = w;

			double x = xyz[0];
			double y = xyz[1];
			double t = xyz[2];

			if (!Double.isFinite(x) || !Double.isFinite(y) || !Double.isFinite(t)) {
				continue;
			}

			double tTol = ROOT_EPS * (Math.abs(t) + 1.0);
			if (t < -tTol) {
				continue;
			}
			if (t < 0.0) {
				t = 0.0;
			}

			if (!satisfiesAll(original, x, y, t)) {
				continue;
			}

			if (isDuplicate(solns, x, y)) {
				continue;
			}

			solns.add(new Solution(new Point(x, y), t, k3));
			added++;
		}

		return added;
	}

	private int solveReducedQuadratic(double u2, double u1, double v2, double v1, double w2, double w1, double g0, double uSlope, double uIntercept,
			double vSlope, double vIntercept, double[] roots) {

		double A = u2 * uSlope * uSlope + v2 * vSlope * vSlope + w2;
		double B = 2.0 * u2 * uSlope * uIntercept + 2.0 * v2 * vSlope * vIntercept + u1 * uSlope + v1 * vSlope + w1;
		double C = u2 * uIntercept * uIntercept + v2 * vIntercept * vIntercept + u1 * uIntercept + v1 * vIntercept + g0;

		return solveQuadratic(A, B, C, roots);
	}

	private int solveQuadratic(double a, double b, double c, double[] roots) {
		double scale = Math.abs(a) + Math.abs(b) + Math.abs(c) + 1.0;

		if (Math.abs(a) <= ROOT_EPS * scale) {
			if (Math.abs(b) <= ROOT_EPS * scale) {
				return 0;
			}
			roots[0] = -c / b;
			return 1;
		}

		double disc = b * b - 4.0 * a * c;
		double discScale = Math.abs(b * b) + Math.abs(4.0 * a * c) + 1.0;

		if (disc < -ROOT_EPS * discScale) {
			return 0;
		}
		if (Math.abs(disc) <= ROOT_EPS * discScale) {
			roots[0] = -b / (2.0 * a);
			return 1;
		}

		double sqrtDisc = Math.sqrt(Math.max(0.0, disc));
		double q = -0.5 * (b + Math.copySign(sqrtDisc, b));

		double r1 = q / a;
		double r2 = (Math.abs(q) <= ROOT_EPS * scale) ? (-b / (2.0 * a)) : (c / q);

		roots[0] = r1;
		if (nearlyEqual(r1, r2, DUP_EPS)) {
			return 1;
		}
		roots[1] = r2;
		return 2;
	}

	private boolean satisfiesAll(List<Eq> eqs, double x, double y, double t) {
		for (var eq : eqs) {
			double residual;
			double scale;

			if (eq.q) {
				residual = x * x + eq.a * x + y * y + eq.b * y - t * t + eq.k * t + eq.c;
				scale = x * x + Math.abs(eq.a * x) + y * y + Math.abs(eq.b * y) + t * t + Math.abs(eq.k * t) + Math.abs(eq.c) + 1.0;
			} else {
				residual = eq.a * x + eq.b * y + eq.k * t + eq.c;
				scale = Math.abs(eq.a * x) + Math.abs(eq.b * y) + Math.abs(eq.k * t) + Math.abs(eq.c) + 1.0;
			}

			if (Math.abs(residual) > RES_EPS * scale) {
				return false;
			}
		}
		return true;
	}

	private boolean isDuplicate(List<Solution> solns, double x, double y) {
		for (var s : solns) {
			if (nearlyEqual(s.p.x, x, DUP_EPS) && nearlyEqual(s.p.y, y, DUP_EPS)) {
				return true;
			}
		}
		return false;
	}

	private int[] permutationOrder(List<Eq> lins) {
		double[] score = new double[PERMUTATIONS.length];
		for (int i = 0; i < PERMUTATIONS.length; i++) {
			int xi = PERMUTATIONS[i][0];
			int yi = PERMUTATIONS[i][1];
			score[i] = permutationScore(lins, xi, yi);
		}

		int[] order = { 0, 1, 2 };
		for (int i = 0; i < order.length; i++) {
			for (int j = i + 1; j < order.length; j++) {
				if (score[order[j]] > score[order[i]]) {
					int tmp = order[i];
					order[i] = order[j];
					order[j] = tmp;
				}
			}
		}
		return order;
	}

	private double permutationScore(List<Eq> lins, int xi, int yi) {
		double ai = lins.get(0).get(xi);
		double bi = lins.get(0).get(yi);
		double aj = lins.get(1).get(xi);
		double bj = lins.get(1).get(yi);
		return Math.abs(ai * bj - aj * bi);
	}

	private boolean nearlyEqual(double a, double b, double eps) {
		double scale = Math.max(1.0, Math.max(Math.abs(a), Math.abs(b)));
		return Math.abs(a - b) <= eps * scale;
	}
}