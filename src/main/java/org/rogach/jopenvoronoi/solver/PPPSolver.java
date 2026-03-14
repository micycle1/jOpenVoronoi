package org.rogach.jopenvoronoi.solver;

import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Point-point-point solver.
 *
 * <p>Computes the circumcenter of the three point sites. The implementation
 * uses a translated 2x2 formulation (relative to the third point) for improved
 * numerical stability. If the points are collinear or nearly collinear, no
 * Voronoi vertex is produced and the solver returns {@code 0}.
 */
public class PPPSolver extends Solver {

	private static final double COLLINEAR_EPS = 1e-12;
	private static final double RES_EPS = 1e-9; // residual tolerance

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		if (!(s1.isPoint() && s2.isPoint() && s3.isPoint())) {
			throw new AssertionError("s1.isPoint() && s2.isPoint() && s3.isPoint()");
		}

		double x1 = s1.x(), y1 = s1.y();
		double x2 = s2.x(), y2 = s2.y();
		double x3 = s3.x(), y3 = s3.y();

		// Translate to p3 to reduce cancellation.
		double ax = x1 - x3;
		double ay = y1 - y3;
		double bx = x2 - x3;
		double by = y2 - y3;

		double a2 = ax * ax + ay * ay;
		double b2 = bx * bx + by * by;

		double det = 2.0 * (ax * by - ay * bx);
		double detScale = 2.0 * (Math.abs(ax * by) + Math.abs(ay * bx)) + 1.0;
		if (!Double.isFinite(det) || Math.abs(det) <= COLLINEAR_EPS * detScale) {
			return 0;
		}

		double uxRel = (by * a2 - ay * b2) / det;
		double uyRel = (ax * b2 - bx * a2) / det;

		double ux = uxRel + x3;
		double uy = uyRel + y3;
		if (!Double.isFinite(ux) || !Double.isFinite(uy)) {
			return 0;
		}

		double dx = ux - x1;
		double dy = uy - y1;
		double t2 = dx * dx + dy * dy;
		if (!(t2 >= 0.0) || !Double.isFinite(t2)) {
			return 0;
		}

		double t = Math.sqrt(t2);
		if (!Double.isFinite(t)) {
			return 0;
		}

		var p = new Point(ux, uy);

		if (!satisfiesPoint(s1, p, t)) {
			return 0;
		}
		if (!satisfiesPoint(s2, p, t)) {
			return 0;
		}
		if (!satisfiesPoint(s3, p, t)) {
			return 0;
		}

		slns.add(new Solution(p, t, +1));
		return 1;
	}

	private boolean satisfiesPoint(Site s, Point p, double t) {
		double dx = p.x - s.x();
		double dy = p.y - s.y();
		double r = dx * dx + dy * dy - t * t;
		double scale = dx * dx + dy * dy + t * t + 1.0;
		return Double.isFinite(r) && Math.abs(r) <= RES_EPS * scale;
	}
}