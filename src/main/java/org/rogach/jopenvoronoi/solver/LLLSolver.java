package org.rogach.jopenvoronoi.solver;

import java.util.ArrayList;
import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Line-line-line solver.
 *
 * <p>
 * Solves the 3x3 linear system
 * 
 * <pre>
 *   a1 x + b1 y + c1 + k1 t = 0
 *   a2 x + b2 y + c2 + k2 t = 0
 *   a3 x + b3 y + c3 + k3 t = 0
 * </pre>
 *
 * <p>
 * The primary solve path uses an allocation-free, unrolled Gaussian elimination
 * with scaled partial pivoting on the unknowns {@code (x, y, t)}. This is more
 * numerically stable than the previous Cramer's-rule-based implementation while
 * avoiding matrix allocation overhead.
 *
 * <p>
 * If the system is singular or nearly singular, the solver falls back to
 * {@link LLLPARASolver} for the special case where two of the three line sites
 * are parallel.
 */
public class LLLSolver extends Solver {

	private static final double PIVOT_EPS = 1e-12;
	private static final double PARALLEL_EPS = 1e-12;
	private static final double T_EPS = 1e-12;
	private static final double RES_EPS = 1e-9;
	private static final double DIR_EPS = 1e-12;

	// output of the allocation-free 3x3 solver
	private double solX;
	private double solY;
	private double solT;

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		if (!(s1.isLine() && s2.isLine() && s3.isLine())) {
			return 0;
		}

		double a1 = s1.a(), b1 = s1.b(), c1 = s1.c();
		double a2 = s2.a(), b2 = s2.b(), c2 = s2.c();
		double a3 = s3.a(), b3 = s3.b(), c3 = s3.c();

		// Solve:
		// [a b k] [x y t]^T = [-c]
		if (solve3x3Unrolled(a1, b1, k1, -c1, a2, b2, k2, -c2, a3, b3, k3, -c3)) {

			double x = solX;
			double y = solY;
			double t = solT;

			if (Double.isFinite(x) && Double.isFinite(y) && Double.isFinite(t)) {
				double tTol = T_EPS * (Math.abs(t) + 1.0);
				if (t >= -tTol) {
					if (t < 0.0) {
						t = 0.0;
					}
					var p = new Point(x, y);
					if (satisfiesAll(a1, b1, c1, k1, a2, b2, c2, k2, a3, b3, c3, k3, p, t)) {
						slns.add(new Solution(p, t, k3));
						return 1;
					}
				}
			}
		}

		return solveParallelFallback(s1, k1, s2, k2, s3, k3, slns);
	}

	/**
	 * Allocation-free, unrolled Gaussian elimination with scaled partial pivoting.
	 *
	 * Solves: [m00 m01 m02] [x] [m03] [m10 m11 m12] [y] = [m13] [m20 m21 m22] [t]
	 * [m23]
	 */
	private boolean solve3x3Unrolled(double m00, double m01, double m02, double m03, double m10, double m11, double m12, double m13, double m20, double m21,
			double m22, double m23) {

		double s0 = rowScale(m00, m01, m02, m03);
		double s1 = rowScale(m10, m11, m12, m13);
		double s2 = rowScale(m20, m21, m22, m23);

		// Pivot column 0
		double p0 = Math.abs(m00) / s0;
		double p1 = Math.abs(m10) / s1;
		double p2 = Math.abs(m20) / s2;

		if (p1 > p0 && p1 >= p2) {
			double t;
			t = m00;
			m00 = m10;
			m10 = t;
			t = m01;
			m01 = m11;
			m11 = t;
			t = m02;
			m02 = m12;
			m12 = t;
			t = m03;
			m03 = m13;
			m13 = t;
			t = s0;
			s0 = s1;
			s1 = t;
		} else if (p2 > p0 && p2 > p1) {
			double t;
			t = m00;
			m00 = m20;
			m20 = t;
			t = m01;
			m01 = m21;
			m21 = t;
			t = m02;
			m02 = m22;
			m22 = t;
			t = m03;
			m03 = m23;
			m23 = t;
			t = s0;
			s0 = s2;
			s2 = t;
		}

		if (Math.abs(m00) <= PIVOT_EPS * s0) {
			return false;
		}

		// Eliminate row 1
		if (m10 != 0.0) {
			double f10 = m10 / m00;
			m11 -= f10 * m01;
			m12 -= f10 * m02;
			m13 -= f10 * m03;
			m10 = 0.0;
		}

		// Eliminate row 2
		if (m20 != 0.0) {
			double f20 = m20 / m00;
			m21 -= f20 * m01;
			m22 -= f20 * m02;
			m23 -= f20 * m03;
			m20 = 0.0;
		}

		// Pivot column 1
		double q1 = Math.abs(m11) / s1;
		double q2 = Math.abs(m21) / s2;
		if (q2 > q1) {
			double t;
			t = m10;
			m10 = m20;
			m20 = t;
			t = m11;
			m11 = m21;
			m21 = t;
			t = m12;
			m12 = m22;
			m22 = t;
			t = m13;
			m13 = m23;
			m23 = t;
			t = s1;
			s1 = s2;
			s2 = t;
		}

		if (Math.abs(m11) <= PIVOT_EPS * s1) {
			return false;
		}

		// Eliminate row 2
		if (m21 != 0.0) {
			double f21 = m21 / m11;
			m22 -= f21 * m12;
			m23 -= f21 * m13;
			m21 = 0.0;
		}

		if (Math.abs(m22) <= PIVOT_EPS * s2) {
			return false;
		}

		// Back-substitution
		double t = m23 / m22;
		double y = (m13 - m12 * t) / m11;
		double x = (m03 - m01 * y - m02 * t) / m00;

		if (!Double.isFinite(x) || !Double.isFinite(y) || !Double.isFinite(t)) {
			return false;
		}

		solX = x;
		solY = y;
		solT = t;
		return true;
	}

	private double rowScale(double a, double b, double c, double d) {
		return Math.abs(a) + Math.abs(b) + Math.abs(c) + Math.abs(d) + 1.0;
	}

	private int solveParallelFallback(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		double d12 = parallelMeasure(s1, s2);
		double d23 = parallelMeasure(s2, s3);
		double d31 = parallelMeasure(s3, s1);

		Site pA, pB, third;
		double kA, kB, kThird;
		double best = d12;

		pA = s1;
		kA = k1;
		pB = s2;
		kB = k2;
		third = s3;
		kThird = k3;

		if (d23 < best) {
			best = d23;
			pA = s2;
			kA = k2;
			pB = s3;
			kB = k3;
			third = s1;
			kThird = k1;
		}
		if (d31 < best) {
			best = d31;
			pA = s3;
			kA = k3;
			pB = s1;
			kB = k1;
			third = s2;
			kThird = k2;
		}

		if (!isNearlyParallel(pA, pB)) {
			return 0;
		}

		var paraSolver = new LLLPARASolver();
		List<Solution> tmp = new ArrayList<>(1);
		int n = paraSolver.solve(pA, kA, pB, kB, third, kThird, tmp);
		if (n <= 0) {
			return 0;
		}

		int added = 0;
		for (var s : tmp) {
			if (!Double.isFinite(s.p.x) || !Double.isFinite(s.p.y) || !Double.isFinite(s.t)) {
				continue;
			}

			double tTol = T_EPS * (Math.abs(s.t) + 1.0);
			if (s.t < -tTol) {
				continue;
			}
			double t = (s.t < 0.0) ? 0.0 : s.t;

			if (!satisfiesLine(s1, k1, s.p, t)) {
				continue;
			}
			if (!satisfiesLine(s2, k2, s.p, t)) {
				continue;
			}
			if (!satisfiesLine(s3, k3, s.p, t)) {
				continue;
			}

			if (!hasProperOffsetDirection(third, kThird, s.p)) {
				continue;
			}

			slns.add(new Solution(s.p, t, k3));
			added++;
		}

		return added;
	}

	private double parallelMeasure(Site sA, Site sB) {
		return Math.abs(sA.a() * sB.b() - sB.a() * sA.b());
	}

	private boolean isNearlyParallel(Site sA, Site sB) {
		double cross = parallelMeasure(sA, sB);
		double scale = Math.abs(sA.a() * sB.b()) + Math.abs(sB.a() * sA.b()) + 1.0;
		return cross <= PARALLEL_EPS * scale;
	}

	private boolean satisfiesAll(double a1, double b1, double c1, double k1, double a2, double b2, double c2, double k2, double a3, double b3, double c3,
			double k3, Point p, double t) {
		return satisfiesLine(a1, b1, c1, k1, p, t) && satisfiesLine(a2, b2, c2, k2, p, t) && satisfiesLine(a3, b3, c3, k3, p, t);
	}

	private boolean satisfiesLine(Site s, double k, Point p, double t) {
		return satisfiesLine(s.a(), s.b(), s.c(), k, p, t);
	}

	private boolean satisfiesLine(double a, double b, double c, double k, Point p, double t) {
		double r = a * p.x + b * p.y + c + k * t;
		double scale = Math.abs(a * p.x) + Math.abs(b * p.y) + Math.abs(c) + Math.abs(k * t) + 1.0;
		return Double.isFinite(r) && Math.abs(r) <= RES_EPS * scale;
	}

	private boolean hasProperOffsetDirection(Site s, double k, Point p) {
		var edge = s.end().sub(s.start());
		var rel = p.sub(s.start());
		double cross = edge.cross(rel);
		double scale = edge.norm() * rel.norm() + 1.0;
		return cross * k >= -DIR_EPS * scale;
	}
}