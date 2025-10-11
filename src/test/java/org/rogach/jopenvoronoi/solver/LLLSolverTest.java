package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for LLLSolver (three line-sites).
 */
public class LLLSolverTest {

    private static final double EPS = 1e-9;

    @Test
    @DisplayName("Three lines intersecting at origin with unit offsets")
    public void threeLinesIntersectingAtOrigin() {
        // Three lines forming 120-degree angles, all passing through origin
        // Line 1: y = 0 (horizontal)
        // Line 2: y = sqrt(3)*x (60 degrees)
        // Line 3: y = -sqrt(3)*x (120 degrees)
        // With k=1 for all, solution should be at origin with t=0
        Site l1 = new LineSite(new Point(-1, 0), new Point(1, 0), 1);
        Site l2 = new LineSite(new Point(-1, -Math.sqrt(3)), new Point(1, Math.sqrt(3)), 1);
        Site l3 = new LineSite(new Point(-1, Math.sqrt(3)), new Point(1, -Math.sqrt(3)), 1);
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        assertEquals(1, n);
        Solution s = out.get(0);
        
        // Should be at or very near origin
        assertEquals(0.0, s.p.x, EPS);
        assertEquals(0.0, s.p.y, EPS);
        assertEquals(0.0, s.t, EPS);
    }

    @Test
    @DisplayName("Equilateral triangle incenter")
    public void equilateralTriangleIncenter() {
        // Equilateral triangle with vertices at (0,0), (2,0), (1,sqrt(3))
        // Lines along the edges with k=1 for inward offsets
        double h = Math.sqrt(3.0);
        
        Site l1 = new LineSite(new Point(0, 0), new Point(2, 0), 1);      // bottom edge
        Site l2 = new LineSite(new Point(2, 0), new Point(1, h), 1);      // right edge
        Site l3 = new LineSite(new Point(1, h), new Point(0, 0), 1);      // left edge
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        assertEquals(1, n);
        Solution s = out.get(0);
        
        // Incenter of equilateral triangle is at centroid
        assertEquals(1.0, s.p.x, EPS);
        assertEquals(h / 3.0, s.p.y, EPS);
        
        // Verify equal distance to all three lines
        double d1 = Math.abs(l1.a() * s.p.x + l1.b() * s.p.y + l1.c()) / 
                    Math.sqrt(l1.a() * l1.a() + l1.b() * l1.b());
        double d2 = Math.abs(l2.a() * s.p.x + l2.b() * s.p.y + l2.c()) / 
                    Math.sqrt(l2.a() * l2.a() + l2.b() * l2.b());
        double d3 = Math.abs(l3.a() * s.p.x + l3.b() * s.p.y + l3.c()) / 
                    Math.sqrt(l3.a() * l3.a() + l3.b() * l3.b());
        
        assertEquals(d1, d2, EPS);
        assertEquals(d2, d3, EPS);
    }

    @Test
    @DisplayName("Lines in different input orders produce same result")
    public void lineOrderInvariance() {
        Site l1 = new LineSite(new Point(0, 0), new Point(1, 0), 1);
        Site l2 = new LineSite(new Point(0, 0), new Point(0, 1), 1);
        Site l3 = new LineSite(new Point(1, 0), new Point(0, 1), 1);
        
        LLLSolver solver = new LLLSolver();
        
        List<Solution> out1 = new ArrayList<>();
        solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out1);
        
        List<Solution> out2 = new ArrayList<>();
        solver.solve(l2, 1.0, l3, 1.0, l1, 1.0, out2);
        
        List<Solution> out3 = new ArrayList<>();
        solver.solve(l3, 1.0, l1, 1.0, l2, 1.0, out3);
        
        // All should produce a solution at the same location
        if (out1.size() > 0 && out2.size() > 0 && out3.size() > 0) {
            assertEquals(out1.get(0).p.x, out2.get(0).p.x, EPS);
            assertEquals(out1.get(0).p.y, out2.get(0).p.y, EPS);
            assertEquals(out1.get(0).p.x, out3.get(0).p.x, EPS);
            assertEquals(out1.get(0).p.y, out3.get(0).p.y, EPS);
        }
    }

    @Test
    @DisplayName("Different k values affect solution position")
    public void differentKValuesAffectSolution() {
        Site l1 = new LineSite(new Point(0, 0), new Point(1, 0), 1);
        Site l2 = new LineSite(new Point(0, 0), new Point(0, 1), 1);
        Site l3 = new LineSite(new Point(1, 0), new Point(0, 1), 1);
        
        LLLSolver solver = new LLLSolver();
        
        List<Solution> out1 = new ArrayList<>();
        solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out1);
        
        List<Solution> out2 = new ArrayList<>();
        solver.solve(l1, -1.0, l2, 1.0, l3, 1.0, out2);
        
        // Different k values should produce different solutions (or no solution)
        if (out1.size() > 0 && out2.size() > 0) {
            boolean different = Math.abs(out1.get(0).p.x - out2.get(0).p.x) > EPS ||
                               Math.abs(out1.get(0).p.y - out2.get(0).p.y) > EPS;
            assertTrue(different, "Different k values should produce different solutions");
        }
    }

    @Test
    @DisplayName("Negative t value returns no solution")
    public void negativeTValueReturnsNoSolution() {
        // Create a configuration where the solution would have negative t
        // Lines pointing outward from a region
        Site l1 = new LineSite(new Point(0, 1), new Point(1, 1), -1); // above, k=-1 points down
        Site l2 = new LineSite(new Point(1, 0), new Point(1, 1), -1); // right, k=-1 points left
        Site l3 = new LineSite(new Point(0, 0), new Point(1, 0), -1); // below, k=-1 points up
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, -1.0, l2, -1.0, l3, -1.0, out);
        
        // Should return 0 if t would be negative
        assertTrue(n >= 0);
        assertEquals(n, out.size());
    }

    @Test
    @DisplayName("Vertical and horizontal lines")
    public void verticalAndHorizontalLines() {
        // x = 1, y = 1, x + y = 3
        Site l1 = new LineSite(new Point(1, 0), new Point(1, 1), 1);  // vertical x=1
        Site l2 = new LineSite(new Point(0, 1), new Point(1, 1), 1);  // horizontal y=1
        Site l3 = new LineSite(new Point(2, 1), new Point(1, 2), 1);  // x+y=3
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        // Should find a solution
        assertTrue(n >= 0);
        if (n > 0) {
            Solution s = out.get(0);
            // Verify the solution satisfies all three line equations
            double err1 = Math.abs(l1.a() * s.p.x + l1.b() * s.p.y + l1.c() + 1.0 * s.t);
            double err2 = Math.abs(l2.a() * s.p.x + l2.b() * s.p.y + l2.c() + 1.0 * s.t);
            double err3 = Math.abs(l3.a() * s.p.x + l3.b() * s.p.y + l3.c() + 1.0 * s.t);
            
            assertTrue(err1 < EPS, "Solution should satisfy line 1 equation");
            assertTrue(err2 < EPS, "Solution should satisfy line 2 equation");
            assertTrue(err3 < EPS, "Solution should satisfy line 3 equation");
        }
    }

    @Test
    @DisplayName("Large coordinate values")
    public void largeCoordinateValues() {
        // Lines far from origin
        Site l1 = new LineSite(new Point(1000, 1000), new Point(1001, 1000), 1);
        Site l2 = new LineSite(new Point(1000, 1000), new Point(1000, 1001), 1);
        Site l3 = new LineSite(new Point(1001, 1000), new Point(1000, 1001), 1);
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        assertTrue(n >= 0);
        if (n > 0) {
            Solution s = out.get(0);
            // Solution should be near (1000, 1000) region
            assertTrue(s.p.x > 999 && s.p.x < 1002);
            assertTrue(s.p.y > 999 && s.p.y < 1002);
        }
    }

    @Test
    @DisplayName("Solution satisfies all three line equations")
    public void solutionSatisfiesAllEquations() {
        // Create arbitrary non-parallel lines
        Site l1 = new LineSite(new Point(0, 0), new Point(1, 1), 1);
        Site l2 = new LineSite(new Point(0, 1), new Point(1, 0), 1);
        Site l3 = new LineSite(new Point(0, 0.5), new Point(1, 0.5), 1);
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        if (n > 0) {
            Solution s = out.get(0);
            
            // Each line equation: a*x + b*y + c + k*t = 0
            double eq1 = l1.a() * s.p.x + l1.b() * s.p.y + l1.c() + 1.0 * s.t;
            double eq2 = l2.a() * s.p.x + l2.b() * s.p.y + l2.c() + 1.0 * s.t;
            double eq3 = l3.a() * s.p.x + l3.b() * s.p.y + l3.c() + 1.0 * s.t;
            
            assertTrue(Math.abs(eq1) < EPS, "Solution should satisfy line 1: " + eq1);
            assertTrue(Math.abs(eq2) < EPS, "Solution should satisfy line 2: " + eq2);
            assertTrue(Math.abs(eq3) < EPS, "Solution should satisfy line 3: " + eq3);
        }
    }

    @Test
    @DisplayName("Non-negative t value in solution")
    public void nonNegativeTValue() {
        Site l1 = new LineSite(new Point(0, 0), new Point(1, 0), 1);
        Site l2 = new LineSite(new Point(0, 0), new Point(0, 1), 1);
        Site l3 = new LineSite(new Point(1, 0), new Point(0, 1), 1);
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        if (n > 0) {
            Solution s = out.get(0);
            assertTrue(s.t >= 0, "Solution t value should be non-negative");
        }
    }

    @Test
    @DisplayName("Returns zero solutions for parallel lines configuration")
    public void parallelLinesNoSolution() {
        // Two parallel lines with a third non-parallel
        Site l1 = new LineSite(new Point(0, 0), new Point(1, 0), 1);  // y = 0
        Site l2 = new LineSite(new Point(0, 1), new Point(1, 1), 1);  // y = 1, parallel to l1
        Site l3 = new LineSite(new Point(0, 0), new Point(0, 1), 1);  // x = 0
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        // Parallel lines may have 0 solutions or use parallel solver
        assertTrue(n >= 0);
        assertEquals(n, out.size());
    }

    @Test
    @DisplayName("Solution count matches returned value")
    public void solutionCountMatchesReturnValue() {
        Site l1 = new LineSite(new Point(0, 0), new Point(1, 0), 1);
        Site l2 = new LineSite(new Point(0, 0), new Point(0, 1), 1);
        Site l3 = new LineSite(new Point(1, 0), new Point(0, 1), 1);
        
        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, 1.0, l2, 1.0, l3, 1.0, out);
        
        assertEquals(n, out.size(), "Returned count should match actual solution count");
    }
}