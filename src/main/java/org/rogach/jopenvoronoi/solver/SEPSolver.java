package org.rogach.jopenvoronoi.solver;

import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Separator solver for the line/point/line case.
 *
 * s1 = line site
 * s2 = endpoint point site of s1
 * s3 = line site
 */
public class SEPSolver extends Solver {

	private static final double DEN_EPS = 1e-12;
	private static final double RES_EPS = 1e-9;

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		if (!(s1.isLine() && s2.isPoint() && s3.isLine())) {
			return 0;
		}

		var sv = separatorDirection(s1, k1);

		double num = -(s3.a() * s2.x() + s3.b() * s2.y() + s3.c());
		double den = sv.x * s3.a() + sv.y * s3.b() + k3;

		double denScale = Math.abs(sv.x * s3.a()) + Math.abs(sv.y * s3.b()) + Math.abs(k3) + 1.0;
		if (!Double.isFinite(den) || Math.abs(den) <= DEN_EPS * denScale) {
			return 0;
		}

		double t = num / den;
		if (!Double.isFinite(t)) {
			return 0;
		}

		double tTol = DEN_EPS * (Math.abs(t) + 1.0);
		if (t < -tTol) {
			return 0;
		}
		if (t < 0.0) {
			t = 0.0;
		}

		var p = new Point(s2.x() + sv.x * t, s2.y() + sv.y * t);

		if (!satisfiesLine(s1, k1, p, t)) {
			return 0;
		}
		if (!satisfiesPoint(s2, p, t)) {
			return 0;
		}
		if (!satisfiesLine(s3, k3, p, t)) {
			return 0;
		}

		slns.add(new Solution(p, t, k3));
		return 1;
	}

	private Point separatorDirection(Site line, double k) {
		// Matches the convention used elsewhere in this codebase:
		// k = +1 -> (-a,-b), k = -1 -> (a,b)
		return (k < 0.0)
				? new Point(line.a(), line.b())
				: new Point(-line.a(), -line.b());
	}

	private boolean satisfiesLine(Site s, double k, Point p, double t) {
		double r = s.a() * p.x + s.b() * p.y + s.c() + k * t;
		double scale = Math.abs(s.a() * p.x) + Math.abs(s.b() * p.y) + Math.abs(s.c()) + Math.abs(k * t) + 1.0;
		return Double.isFinite(r) && Math.abs(r) <= RES_EPS * scale;
	}

	private boolean satisfiesPoint(Site s, Point p, double t) {
		double dx = p.x - s.x();
		double dy = p.y - s.y();
		double r = dx * dx + dy * dy - t * t;
		double scale = dx * dx + dy * dy + t * t + 1.0;
		return Double.isFinite(r) && Math.abs(r) <= RES_EPS * scale;
	}
}