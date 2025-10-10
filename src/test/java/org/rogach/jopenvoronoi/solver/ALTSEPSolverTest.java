package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;
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
 * Tests for ALTSEPSolver (alternative separator cases).
 */
public class ALTSEPSolverTest {

	@Test
	public void altSeparatorType0MatchesExpected() {
		// s3 (line) and s1 (point) form the separator pair for type=0.
		Site p = new PointSite(new Point(0, 0)); // s1 (point)
		Site lineX1 = new LineSite(new Point(1, -1), new Point(1, 1), +1); // s2 (third site): x=1
		Site lineY0 = new LineSite(new Point(0, 0), new Point(1, 0), -1); // s3 (separator line): y=0, k=-1

		ALTSEPSolver solver = new ALTSEPSolver();
		solver.set_type(0);

		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p, 1, lineX1, lineX1.k(), lineY0, lineY0.k(), out);

		assertEquals(1, n);
		Solution sol = out.get(0);

		// Position on separator x=0, at unit distance from the point (0,0)
		assertEquals(0.0, sol.p.x, 1e-12);
		assertEquals(1.0, Math.abs(sol.p.y), 1e-12);

		// Offset distance magnitude should be 1
		assertEquals(1.0, Math.abs(sol.t), 1e-12);

		// verify third-site constraint a3*x + b3*y + c3 + t ≈ 0
		double a3 = lineX1.a(), b3 = lineX1.b(), c3 = lineX1.c();
		double lhs = a3 * sol.p.x + b3 * sol.p.y + c3 + sol.t;
		assertEquals(0.0, lhs, 1e-12);
	}

	@Test
	public void altSeparatorType0WithPointThirdSite() {
		Site p = new PointSite(new Point(0, 0)); // s1 (point)
		Site thirdPoint = new PointSite(new Point(0, 1)); // s2 (third site: point)
		Site lineY0 = new LineSite(new Point(0, 0), new Point(1, 0), -1); // s3 (separator line): y=0, k=-1

		ALTSEPSolver solver = new ALTSEPSolver();
		solver.set_type(0);

		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p, 1, thirdPoint, 1, lineY0, lineY0.k(), out);

		assertEquals(1, n);
		Solution sol = out.get(0);

		assertEquals(0.0, sol.p.x, 1e-12);
		assertEquals(0.5, sol.p.y, 1e-12);
		assertEquals(0.5, Math.abs(sol.t), 1e-12);
		assertEquals(lineY0.k(), sol.k3, 0.0);

		double distToThird = Math.hypot(sol.p.x - thirdPoint.x(), sol.p.y - thirdPoint.y());
		assertEquals(Math.abs(sol.t), distToThird, 1e-12);
	}

	@Test
	public void altSeparatorType1WithLineThirdSite() {
		Site p = new PointSite(new Point(0, 0)); // s2 (point)
		Site lineX1 = new LineSite(new Point(1, -1), new Point(1, 1), +1); // s1 (third site): x=1
		Site lineY0 = new LineSite(new Point(0, 0), new Point(1, 0), -1); // s3 (separator line): y=0, k=-1

		ALTSEPSolver solver = new ALTSEPSolver();
		solver.set_type(1);

		List<Solution> out = new ArrayList<>();
		int n = solver.solve(lineX1, lineX1.k(), p, 1, lineY0, lineY0.k(), out);

		assertEquals(1, n);
		Solution sol = out.get(0);

		assertEquals(0.0, sol.p.x, 1e-12);
		assertEquals(1.0, Math.abs(sol.p.y), 1e-12);
		assertEquals(1.0, Math.abs(sol.t), 1e-12);
		assertEquals(lineY0.k(), sol.k3, 0.0);

		// Check the line constraint uses +t as in solver
		double lhs = lineX1.a() * sol.p.x + lineX1.b() * sol.p.y + lineX1.c() + sol.t;
		assertEquals(0.0, lhs, 1e-12);
	}

	@Test
	public void altSeparatorType0NoSolutionPointCase() {
		Site p = new PointSite(new Point(0, 0)); // s1 (point)
		Site thirdPoint = new PointSite(new Point(1, 0)); // s2 (third site)
		Site lineY0 = new LineSite(new Point(0, 0), new Point(1, 0), -1); // s3: y=0, k=-1 -> sv=(0,1)

		ALTSEPSolver solver = new ALTSEPSolver();
		solver.set_type(0);

		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p, 1, thirdPoint, 1, lineY0, lineY0.k(), out);

		assertEquals(0, n);
		assertTrue(out.isEmpty());
	}

	@Test
	public void altSeparatorType0NoSolutionLineCase() {
		Site p = new PointSite(new Point(0, 0)); // s1
		Site lineY0 = new LineSite(new Point(0, 0), new Point(1, 0), -1); // s3: y=0, k=-1 -> sv=(0,1)
		// third line with b ≈ -1: reverse endpoints of a horizontal line
		Site thirdLine = new LineSite(new Point(1, 2), new Point(0, 2), +1);

		ALTSEPSolver solver = new ALTSEPSolver();
		solver.set_type(0);

		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p, 1, thirdLine, thirdLine.k(), lineY0, lineY0.k(), out);

		// If thirdLine.b() == -1, then n == 0; otherwise, you can assert near-zero
		// denominator and skip.
		if (Math.abs(lineY0.a() * 0 + 1 * thirdLine.b() + 1) < 1e-12) {
			assertEquals(0, n);
			assertTrue(out.isEmpty());
		}
	}

}
