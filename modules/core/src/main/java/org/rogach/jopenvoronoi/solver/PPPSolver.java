package org.rogach.jopenvoronoi.solver;

import static org.rogach.jopenvoronoi.util.Numeric.sq;

import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * point-point-point Solver (based on Sugihara &amp; Iri paper)
 */
public class PPPSolver extends Solver {

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {
		assert (s1.isPoint() && s2.isPoint() && s3.isPoint()) : "s1.isPoint() && s2.isPoint() && s3.isPoint()";
		var pi = s1.position();
		var pj = s2.position();
		var pk = s3.position();

		if (pi.isRight(pj, pk)) {
			var tmp = pi;
			pi = pj;
			pj = tmp;
		}
		assert (!pi.isRight(pj, pk)) : " !pi.is_right(pj,pk) ";
		// 2) point pk should have the largest angle. largest angle is opposite longest
		// side.
		// squared distances: sqrt is monotonic, so comparisons are identical to
		// comparing distance(); no value from this loop flows into the output
		var longest_side_sq = distanceSq(pi, pj);
		while ((distanceSq(pj, pk) > longest_side_sq) || ((distanceSq(pi, pk) > longest_side_sq))) {
			// cyclic rotation of points until pk is opposite the longest side pi-pj
			var tmp = pk;
			pk = pj;
			pj = pi;
			pi = tmp;
			longest_side_sq = distanceSq(pi, pj);
		}
		assert (!pi.isRight(pj, pk)) : " !pi.is_right(pj,pk) ";
		assert (pi.sub(pj).norm() >= pj.sub(pk).norm()) : " pi.sub(pj).norm() >=  pj.sub(pk).norm() ";
		assert (pi.sub(pj).norm() >= pk.sub(pi).norm()) : " pi.sub(pj).norm() >=  pk.sub(pi).norm() ";

		var J2 = (pi.y - pk.y) * (sq(pj.x - pk.x) + sq(pj.y - pk.y)) / 2.0
				- (pj.y - pk.y) * (sq(pi.x - pk.x) + sq(pi.y - pk.y)) / 2.0;
		var J3 = (pi.x - pk.x) * (sq(pj.x - pk.x) + sq(pj.y - pk.y)) / 2.0
				- (pj.x - pk.x) * (sq(pi.x - pk.x) + sq(pi.y - pk.y)) / 2.0;
		var J4 = (pi.x - pk.x) * (pj.y - pk.y) - (pj.x - pk.x) * (pi.y - pk.y);
		assert (J4 != 0.0) : " J4 != 0.0 ";
		if (J4 == 0.0) {
			throw new RuntimeException(" PPPSolver: Warning divide-by-zero!!");
		}
		var sln_pt = new Point(-J2 / J4 + pk.x, J3 / J4 + pk.y);
		var dist = sln_pt.distance(pi);
		slns.add(new Solution(sln_pt, dist, +1));
		return 1;
	}

	/** squared distance; exactly the argument of the sqrt in {@link Point#distance(Point)} */
	private static double distanceSq(Point a, Point b) {
		var dx = a.x - b.x;
		var dy = a.y - b.y;
		return dx * dx + dy * dy;
	}

}
