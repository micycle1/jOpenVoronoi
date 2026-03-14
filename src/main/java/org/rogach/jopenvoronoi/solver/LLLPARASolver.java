package org.rogach.jopenvoronoi.solver;

import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.util.Pair;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Line-line-line solver for the parallel line-segment case.
 *
 * <p>
 * Used when two of the three line sites are parallel. The solution is obtained
 * by intersecting:
 * <ul>
 * <li>the line-line bisector of the parallel pair, and</li>
 * <li>the offset line of the third site at distance {@code t_b}.</li>
 * </ul>
 *
 * <p>
 * The 2x2 system is solved with a scaled determinant test for improved
 * robustness.
 */
public class LLLPARASolver extends Solver {

	private static final double PARALLEL_EPS = 1e-12;
	private static final double DET_EPS = 1e-12;
	private static final double RES_EPS = 1e-9;

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {

		if (!isNearlyParallel(s1, s2)) {
			return 0;
		}

		double ba = s1.a();
		double bb = s1.b();
		double s2c = s2.c();

		// If normals point in opposite directions, align the second equation first.
		double dot = s1.a() * s2.a() + s1.b() * s2.b();
		if (dot < 0.0) {
			s2c = -s2c;
		}

		double bc = 0.5 * (s1.c() + s2c);
		double tb = 0.5 * Math.abs(s1.c() - s2c);

		var xy = twoByTwoSolve(ba, bb, s3.a(), s3.b(), -bc, -s3.c() - k3 * tb);

		if (xy == null) {
			return 0;
		}

		double x = xy.getFirst();
		double y = xy.getSecond();
		if (!Double.isFinite(x) || !Double.isFinite(y) || !Double.isFinite(tb)) {
			return 0;
		}

		var p = new Point(x, y);

		// Validate against the original three offset equations.
		if (!satisfiesLine(s1, k1, p, tb)) {
			return 0;
		}
		if (!satisfiesLine(s2, k2, p, tb)) {
			return 0;
		}
		if (!satisfiesLine(s3, k3, p, tb)) {
			return 0;
		}

		slns.add(new Solution(p, tb, k3));
		return 1;
	}

	private boolean isNearlyParallel(Site s1, Site s2) {
		double cross = Math.abs(s1.a() * s2.b() - s2.a() * s1.b());
		double scale = Math.abs(s1.a() * s2.b()) + Math.abs(s2.a() * s1.b()) + 1.0;
		return cross <= PARALLEL_EPS * scale;
	}

	private Pair<Double, Double> twoByTwoSolve(double a, double b, double c, double d, double e, double f) {
		// [ a b ] [u] = [ e ]
		// [ c d ] [v] = [ f ]
		double det = a * d - c * b;
		double scale = Math.abs(a * d) + Math.abs(c * b) + 1.0;
		if (!Double.isFinite(det) || Math.abs(det) <= DET_EPS * scale) {
			return null;
		}

		double u = (d * e - b * f) / det;
		double v = (-c * e + a * f) / det;
		if (!Double.isFinite(u) || !Double.isFinite(v)) {
			return null;
		}

		return new Pair<>(u, v);
	}

	private boolean satisfiesLine(Site s, double k, Point p, double t) {
		double r = s.a() * p.x + s.b() * p.y + s.c() + k * t;
		double scale = Math.abs(s.a() * p.x) + Math.abs(s.b() * p.y) + Math.abs(s.c()) + Math.abs(k * t) + 1.0;
		return Double.isFinite(r) && Math.abs(r) <= RES_EPS * scale;
	}
}