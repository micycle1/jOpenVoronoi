package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

public class PPPSolverTest {

	private static final double EPS = 1e-9;

	@Test
	@DisplayName("Circumcenter of right triangle")
	public void circumcenterOfRightTriangle() {
		// Right triangle points (0,0), (1,0), (0,1)
		// Circumcenter at (0.5, 0.5) with radius sqrt(0.5)
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(1, 0));
		Site p3 = new PointSite(new Point(0, 1));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);
		assertEquals(0.5, s.p.x, EPS);
		assertEquals(0.5, s.p.y, EPS);
		assertEquals(Math.sqrt(0.5), s.t, EPS);
	}

	@Test
	@DisplayName("Equilateral triangle - circumcenter at centroid")
	public void equilateralTriangle() {
		// Equilateral triangle with vertices at specific positions
		double h = Math.sqrt(3.0) / 2.0; // height of unit equilateral triangle
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(1, 0));
		Site p3 = new PointSite(new Point(0.5, h));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);

		// Circumcenter should be at centroid (1/3, h/3) with radius 1/sqrt(3)
		double expectedX = 0.5;
		double expectedY = h / 3.0;
		double expectedRadius = 1.0 / Math.sqrt(3.0);

		assertEquals(expectedX, s.p.x, EPS);
		assertEquals(expectedY, s.p.y, EPS);
		assertEquals(expectedRadius, s.t, EPS);
	}

	@Test
	@DisplayName("Isosceles triangle")
	public void isoscelesTriangle() {
		// Isosceles triangle: (0,0), (2,0), (1,2)
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(2, 0));
		Site p3 = new PointSite(new Point(1, 2));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);

		// Circumcenter should be on the axis of symmetry (x=1)
		assertEquals(1.0, s.p.x, EPS);

		// Verify equal distances to all three points
		double d1 = s.p.sub(p1.position()).norm();
		double d2 = s.p.sub(p2.position()).norm();
		double d3 = s.p.sub(p3.position()).norm();

		assertEquals(d1, d2, EPS);
		assertEquals(d2, d3, EPS);
		assertEquals(d1, s.t, EPS);
	}

	@Test
	@DisplayName("Points in different input orders produce same result")
	public void pointOrderInvariance() {
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(4, 0));
		Site p3 = new PointSite(new Point(0, 3));

		PPPSolver solver = new PPPSolver();

		// Test all 6 permutations
		List<Solution> out1 = new ArrayList<>();
		solver.solve(p1, 1, p2, 1, p3, 1, out1);

		List<Solution> out2 = new ArrayList<>();
		solver.solve(p1, 1, p3, 1, p2, 1, out2);

		List<Solution> out3 = new ArrayList<>();
		solver.solve(p2, 1, p1, 1, p3, 1, out3);

		List<Solution> out4 = new ArrayList<>();
		solver.solve(p2, 1, p3, 1, p1, 1, out4);

		List<Solution> out5 = new ArrayList<>();
		solver.solve(p3, 1, p1, 1, p2, 1, out5);

		List<Solution> out6 = new ArrayList<>();
		solver.solve(p3, 1, p2, 1, p1, 1, out6);

		// All should produce the same circumcenter
		assertEquals(out1.get(0).p.x, out2.get(0).p.x, EPS);
		assertEquals(out1.get(0).p.y, out2.get(0).p.y, EPS);
		assertEquals(out1.get(0).p.x, out3.get(0).p.x, EPS);
		assertEquals(out1.get(0).p.y, out3.get(0).p.y, EPS);
		assertEquals(out1.get(0).p.x, out4.get(0).p.x, EPS);
		assertEquals(out1.get(0).p.y, out4.get(0).p.y, EPS);
		assertEquals(out1.get(0).p.x, out5.get(0).p.x, EPS);
		assertEquals(out1.get(0).p.y, out5.get(0).p.y, EPS);
		assertEquals(out1.get(0).p.x, out6.get(0).p.x, EPS);
		assertEquals(out1.get(0).p.y, out6.get(0).p.y, EPS);
	}

	@Test
	@DisplayName("Obtuse triangle")
	public void obtuseTriangle() {
		// Triangle with obtuse angle at p1
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(1, 0));
		Site p3 = new PointSite(new Point(0.2, 0.2));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);

		// Verify equal distances to all three points
		double d1 = s.p.sub(p1.position()).norm();
		double d2 = s.p.sub(p2.position()).norm();
		double d3 = s.p.sub(p3.position()).norm();

		assertEquals(d1, d2, EPS);
		assertEquals(d2, d3, EPS);
		assertEquals(d1, s.t, EPS);
	}

	@Test
	@DisplayName("Large coordinate values")
	public void largeCoordinates() {
		Site p1 = new PointSite(new Point(1000, 1000));
		Site p2 = new PointSite(new Point(1003, 1000));
		Site p3 = new PointSite(new Point(1000, 1004));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);

		// Verify equal distances to all three points
		double d1 = s.p.sub(p1.position()).norm();
		double d2 = s.p.sub(p2.position()).norm();
		double d3 = s.p.sub(p3.position()).norm();

		assertEquals(d1, d2, EPS);
		assertEquals(d2, d3, EPS);
	}

	@Test
	@DisplayName("Negative coordinates")
	public void negativeCoordinates() {
		Site p1 = new PointSite(new Point(-2, -2));
		Site p2 = new PointSite(new Point(-1, -2));
		Site p3 = new PointSite(new Point(-2, -1));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);

		assertEquals(-1.5, s.p.x, EPS);
		assertEquals(-1.5, s.p.y, EPS);
	}

	@Test
	@DisplayName("Very small triangle")
	public void verySmallTriangle() {
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(0.001, 0));
		Site p3 = new PointSite(new Point(0, 0.001));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		Solution s = out.get(0);

		// Verify equal distances to all three points
		double d1 = s.p.sub(p1.position()).norm();
		double d2 = s.p.sub(p2.position()).norm();
		double d3 = s.p.sub(p3.position()).norm();

		assertEquals(d1, d2, EPS);
		assertEquals(d2, d3, EPS);
	}

	@Test
	@DisplayName("Solution offset equals circumradius")
	public void solutionOffsetEqualsCircumradius() {
		Site p1 = new PointSite(new Point(1, 1));
		Site p2 = new PointSite(new Point(4, 1));
		Site p3 = new PointSite(new Point(1, 5));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		solver.solve(p1, 1, p2, 1, p3, 1, out);

		Solution s = out.get(0);

		// s.t should equal the distance from s.p to any of the input points
		double d1 = s.p.sub(p1.position()).norm();
		assertEquals(s.t, d1, EPS);
	}

	@Test
	@DisplayName("Returns exactly one solution")
	public void returnsOneSolution() {
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(5, 0));
		Site p3 = new PointSite(new Point(2, 3));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		int n = solver.solve(p1, 1, p2, 1, p3, 1, out);

		assertEquals(1, n);
		assertEquals(1, out.size());
	}

	@Test
	@DisplayName("Solution has correct sign (+1)")
	public void solutionHasCorrectSign() {
		Site p1 = new PointSite(new Point(0, 0));
		Site p2 = new PointSite(new Point(1, 0));
		Site p3 = new PointSite(new Point(0, 1));

		PPPSolver solver = new PPPSolver();
		List<Solution> out = new ArrayList<>();
		solver.solve(p1, 1, p2, 1, p3, 1, out);

		Solution s = out.get(0);
		assertEquals(1, s.k3);
	}
}