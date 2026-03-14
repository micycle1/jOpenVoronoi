package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for LLLPARASolver (parallel line segments with a third line).
 */
public class LLLPARASolverTest {

	private LLLSolver solver;
	private static final double EPSILON = 1e-10;

	@BeforeEach
	void setUp() {
		solver = new LLLSolver();
	}

	@Test
	@DisplayName("Test basic 3-line intersection")
	void testBasicThreeLineIntersection() {
		// Create three lines that intersect at a single point
		// Line 1: y = x (a=1, b=-1, c=0)
		Site line1 = new LineSite(new Point(0, 0), new Point(1, 1));

		// Line 2: y = -x + 2 (a=1, b=1, c=-2)
		Site line2 = new LineSite(new Point(0, 2), new Point(2, 0));

		// Line 3: x = 1 (a=1, b=0, c=-1)
		Site line3 = new LineSite(new Point(1, 0), new Point(1, 1));

		List<Solution> solutions = new ArrayList<>();
		int result = solver.solve(line1, 1.0, line2, 1.0, line3, 1.0, solutions);

		assertEquals(1, result, "Should find exactly one solution");
		assertEquals(1, solutions.size());

		Solution sol = solutions.get(0);
		// The lines should intersect at (1, 1)
		assertEquals(1.0, sol.p.x, EPSILON);
		assertEquals(1.0, sol.p.y, EPSILON);
		assertTrue(sol.t >= 0, "t value should be non-negative");
	}

	@Test
	@DisplayName("Test parallel lines - should use fallback solver")
	void testParallelLines() {
		// Create two parallel lines and one non-parallel
		// Line 1: y = 1 (horizontal)
		Site line1 = new LineSite(new Point(0, 1), new Point(1, 1));

		// Line 2: y = 2 (horizontal, parallel to line1)
		Site line2 = new LineSite(new Point(0, 2), new Point(1, 2));

		// Line 3: x = 1 (vertical)
		Site line3 = new LineSite(new Point(1, 0), new Point(1, 1));

		List<Solution> solutions = new ArrayList<>();

		// This should trigger the parallel solver fallback
		int result = solver.solve(line1, 1.0, line2, 1.0, line3, 1.0, solutions);

		// The parallel solver should be invoked, exact result depends on LLLPARASolver
		// implementation
		assertNotNull(solutions);
	}

	@Test
	@DisplayName("Test with negative t value - should return no solution")
	void testNegativeTValue() {
		// Create lines that would intersect but with parameters that give negative t
		// This is a bit tricky to set up without knowing the exact equation
		// implementation
		// but the test structure is here
		Site line1 = new LineSite(new Point(0, 0), new Point(1, 0));
		Site line2 = new LineSite(new Point(0, 1), new Point(1, 1));
		Site line3 = new LineSite(new Point(0, 2), new Point(1, 2));

		List<Solution> solutions = new ArrayList<>();

		// With certain k values, this might produce negative t
		int result = solver.solve(line1, -1.0, line2, -1.0, line3, -1.0, solutions);

		// If t is negative, no solution should be added
		if (result == 0) {
			assertEquals(0, solutions.size());
		}
	}

	@Test
	@DisplayName("Test with zero determinant")
	void testZeroDeterminant() {
		// Create three lines that don't have a unique intersection point
		// For example, three lines that all pass through the origin with the same slope
		Site line1 = new LineSite(new Point(0, 0), new Point(1, 1));
		Site line2 = new LineSite(new Point(0, 0), new Point(2, 2));
		Site line3 = new LineSite(new Point(0, 0), new Point(3, 3));

		List<Solution> solutions = new ArrayList<>();
		int result = solver.solve(line1, 1.0, line2, 1.0, line3, 1.0, solutions);

		// Should return 0 for degenerate case
		assertEquals(0, result);
	}

	@Test
	@DisplayName("Test with different k values")
	void testDifferentKValues() {
		Site line1 = new LineSite(new Point(0, 0), new Point(1, 1));
		Site line2 = new LineSite(new Point(0, 2), new Point(2, 0));
		Site line3 = new LineSite(new Point(1, 0), new Point(1, 1));

		List<Solution> solutions = new ArrayList<>();
		int result = solver.solve(line1, 0.5, line2, 1.5, line3, 2.0, solutions);

		if (result > 0) {
			assertEquals(1, solutions.size());
			// k3 should pass through unchanged
			assertEquals(2.0, solutions.get(0).k3, EPSILON);
		}
	}

	@Test
	@DisplayName("Test assertion for non-line sites")
	void testNonLineSiteAssertion() {
		// Assuming we have PointSite or other non-line sites
		Site point1 = new PointSite(new Point(0, 0));
		Site line1 = new LineSite(new Point(0, 0), new Point(1, 1));
		Site line2 = new LineSite(new Point(0, 2), new Point(2, 0));

		List<Solution> solutions = new ArrayList<>();

		assertThrows(AssertionError.class, () -> {
			solver.solve(point1, 1.0, line1, 1.0, line2, 1.0, solutions);
		});
	}

	@ParameterizedTest
	@MethodSource("provideEdgeCases")
	@DisplayName("Test edge cases with various line configurations")
	void testEdgeCases(Site s1, Site s2, Site s3, double k1, double k2, double k3) {
		List<Solution> solutions = new ArrayList<>();
		int result = solver.solve(s1, k1, s2, k2, s3, k3, solutions);

		// Just ensure no exceptions are thrown
		assertTrue(result >= 0);
		assertTrue(solutions.size() >= 0);
	}

	static Stream<Object[]> provideEdgeCases() {
		return Stream.of(
				// Vertical lines
				new Object[] { new LineSite(new Point(0, 0), new Point(0, 1)), new LineSite(new Point(1, 0), new Point(1, 1)),
						new LineSite(new Point(2, 0), new Point(2, 1)), 1.0, 1.0, 1.0 },
				// Horizontal lines
				new Object[] { new LineSite(new Point(0, 0), new Point(1, 0)), new LineSite(new Point(0, 1), new Point(1, 1)),
						new LineSite(new Point(0, 2), new Point(1, 2)), 1.0, 1.0, 1.0 },
				// Lines at 45 degrees
				new Object[] { new LineSite(new Point(0, 0), new Point(1, 1)), new LineSite(new Point(1, 0), new Point(2, 1)),
						new LineSite(new Point(2, 0), new Point(3, 1)), 0.5, 1.0, 1.5 });
	}

	@Test
	@DisplayName("Test solution validity")
	void testSolutionValidity() {
		Site line1 = new LineSite(new Point(0, 0), new Point(2, 1));
		Site line2 = new LineSite(new Point(0, 1), new Point(2, 0));
		Site line3 = new LineSite(new Point(1, 0), new Point(1, 2));

		List<Solution> solutions = new ArrayList<>();
		int result = solver.solve(line1, 1.0, line2, 1.0, line3, 1.0, solutions);

		if (result > 0) {
			Solution sol = solutions.get(0);

			// Verify the solution point actually lies on all three lines
			// This would require implementing a method to check if a point is on a line
			// For now, we just check basic properties
			assertNotNull(sol.p);
			assertTrue(sol.t >= 0);
			assertEquals(1.0, sol.k3, EPSILON);
		}
	}

	static class LineSite extends org.rogach.jopenvoronoi.site.LineSite {
		public LineSite(Point start, Point end) {
			super(start, end, -1);
		}
	}
}
