package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for QLLSolver (quadratic-linear-linear general solver).
 */
public class QLLSolverTest {

	@Test
	public void twoLinesOnePoint() {
		// Lines x=0 and y=0 as segments through the origin
		Site lx0 = new LineSite(new Point(0, 0), new Point(0, 1), -1); // vertical through x=0
		Site ly0 = new LineSite(new Point(0, 0), new Point(1, 0), +1); // horizontal through y=0
		Site p = new PointSite(new Point(1, 1));

		QLLSolver solver = new QLLSolver();

		// Try all k-sign combinations for the two lines
		double[] ks = new double[] { -1.0, +1.0 };
		List<Solution> out = new ArrayList<>();

		for (double kx : ks) {
			for (double ky : ks) {
				solver.solve(lx0, kx, ly0, ky, p, 1.0, out);
			}
		}

		// Tolerances
		double tolLine = 1e-6;
		double tolCirc = 1e-6;

		// Pull coefficients for the line equations
		double a1 = lx0.a(), b1 = lx0.b(), c1 = lx0.c();
		double a2 = ly0.a(), b2 = ly0.b(), c2 = ly0.c();

		// Check if any produced solution satisfies the constraints for some k-pair we
		// tried
		boolean found = false;
		for (Solution s : out) {
			double x = s.p.x;
			double y = s.p.y;
			double t = s.t;

			// See if this solution satisfies the constraints for any of the k-pairs we used
			for (double kx : ks) {
				for (double ky : ks) {
					double r1 = a1 * x + b1 * y + c1 + kx * t;
					double r2 = a2 * x + b2 * y + c2 + ky * t;
					double r3 = (x - 1.0) * (x - 1.0) + (y - 1.0) * (y - 1.0) - t * t;

					if (Math.abs(r1) < tolLine && Math.abs(r2) < tolLine && Math.abs(r3) < tolCirc) {
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
			if (found)
				break;
		}

		assertTrue(found, "Expected at least one QLL solution that satisfies both line equations and the point circle.");

		// Optional: If normals are unit (typical in this codebase), the |t| magnitudes
		// are 2 ± sqrt(2).
		// We won’t fail the test on this; it’s just an extra sanity check.
		double norm1 = Math.hypot(a1, b1);
		double norm2 = Math.hypot(a2, b2);
		if (Math.abs(norm1 - 1.0) < 1e-9 && Math.abs(norm2 - 1.0) < 1e-9) {
			double tSmall = 2.0 - Math.sqrt(2.0);
			double tLarge = 2.0 + Math.sqrt(2.0);
			boolean hasExpectedT = false;
			for (Solution s : out) {
				double tAbs = Math.abs(s.t);
				if (Math.abs(tAbs - tSmall) < 1e-5 || Math.abs(tAbs - tLarge) < 1e-5) {
					hasExpectedT = true;
					break;
				}
			}
			assertTrue(hasExpectedT, "Expected |t| to match one of {2 ± sqrt(2)} when normals are unit length.");
		}
	}
}