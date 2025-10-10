package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for LLLSolver (three line-sites).
 */
public class LLLSolverTest {

    @Test
    public void solvesIncenterForRightTriangleLines() {
        // Lines: x = 0, y = 0, x + y = 1
        // Choose k values so that t equals the (signed) distance to each line and is positive.
        // For x=0 line, use k=-1 so: x - t = 0 -> t = x.
        // For y=0 line, use k=+1 so: -y + t = 0 -> t = y.
        // For x+y=1, use k=+1 so: (x+y-1)/sqrt(2) + t = 0 -> t = (1-(x+y))/sqrt(2).
        // The common solution is the incenter with r = 1/(2+sqrt(2)) at (r, r).

        Site l1 = new LineSite(new Point(0, 0), new Point(0, 1), -1); // x = 0, k = -1
        Site l2 = new LineSite(new Point(0, 0), new Point(1, 0), +1); // y = 0, k = +1
        Site l3 = new LineSite(new Point(1, 0), new Point(0, 1), +1); // x + y = 1, k = +1

        LLLSolver solver = new LLLSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(l1, l1.k(), l2, l2.k(), l3, l3.k(), out);

        assertEquals(1, n, "Expect exactly one solution for three non-parallel lines");
        Solution s = out.get(0);

        double r = 1.0 / (2.0 + Math.sqrt(2.0));
        double eps = 1e-9;
        assertEquals(r, s.p.x, 1e-8);
        assertEquals(r, s.p.y, 1e-8);
        assertEquals(r, s.t, 1e-8);

        // Verify it satisfies each line's offset equation approximately
        assertTrue(Math.abs(l1.a() * s.p.x + l1.b() * s.p.y + l1.c() + l1.k() * s.t) < eps);
        assertTrue(Math.abs(l2.a() * s.p.x + l2.b() * s.p.y + l2.c() + l2.k() * s.t) < eps);
        assertTrue(Math.abs(l3.a() * s.p.x + l3.b() * s.p.y + l3.c() + l3.k() * s.t) < eps);
    }
}

