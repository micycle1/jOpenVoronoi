package org.rogach.jopenvoronoi.solver;

import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Alternative separator solver.
 *
 * type = 0 : s3 / s1 form the separator, s2 is the third site
 * type = 1 : s3 / s2 form the separator, s1 is the third site
 */
public class ALTSEPSolver extends Solver {

	private static final double DEN_EPS = 1e-12;
	private static final double RES_EPS = 1e-9;

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		Site lsite;
		Site psite;
		Site thirdSite;
		double lsiteK;
		double thirdSiteK;

		if (type == 0) {
			lsite = s3;
			lsiteK = k3;
			psite = s1;
			thirdSite = s2;
			thirdSiteK = k2;
		} else if (type == 1) {
			lsite = s3;
			lsiteK = k3;
			psite = s2;
			thirdSite = s1;
			thirdSiteK = k1;
		} else {
			return 0;
		}

		if (!(lsite.isLine() && psite.isPoint())) {
			return 0;
		}
		if (!(thirdSite.isPoint() || thirdSite.isLine())) {
			return 0;
		}

		var sv = separatorDirection(lsite, lsiteK);

		double t;
		if (thirdSite.isPoint()) {
			double dx = psite.x() - thirdSite.x();
			double dy = psite.y() - thirdSite.y();
			double den = 2.0 * (dx * sv.x + dy * sv.y);
			double denScale = 2.0 * (Math.abs(dx * sv.x) + Math.abs(dy * sv.y)) + 1.0;
			if (!Double.isFinite(den) || Math.abs(den) <= DEN_EPS * denScale) {
				return 0;
			}
			t = -(dx * dx + dy * dy) / den;
		} else {
			double num = -(thirdSite.a() * psite.x() + thirdSite.b() * psite.y() + thirdSite.c());
			double den = sv.x * thirdSite.a() + sv.y * thirdSite.b() + thirdSiteK;
			double denScale = Math.abs(sv.x * thirdSite.a()) + Math.abs(sv.y * thirdSite.b()) + Math.abs(thirdSiteK) + 1.0;
			if (!Double.isFinite(den) || Math.abs(den) <= DEN_EPS * denScale) {
				return 0;
			}
			t = num / den;
		}

		if (!Double.isFinite(t)) {
			return 0;
		}

		var p = new Point(psite.x() + sv.x * t, psite.y() + sv.y * t);

		if (!satisfiesLine(lsite, lsiteK, p, t)) {
			return 0;
		}
		if (!satisfiesPoint(psite, p, t)) {
			return 0;
		}
		if (thirdSite.isPoint()) {
			if (!satisfiesPoint(thirdSite, p, t)) {
				return 0;
			}
		} else {
			if (!satisfiesLine(thirdSite, thirdSiteK, p, t)) {
				return 0;
			}
		}

		slns.add(new Solution(p, t, lsiteK));
		return 1;
	}

	private Point separatorDirection(Site line, double k) {
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