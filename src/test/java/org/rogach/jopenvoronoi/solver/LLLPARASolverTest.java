package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for LLLPARASolver (parallel line segments with a third line).
 */
public class LLLPARASolverTest {

    @Test
    public void parallelLinesWithThirdLine() {
        // Parallel lines: y=0 and y=1; third line: x=0 (k3 = -1)
        // Bisector is y=0.5 and tb = 0.5, solve intersection with x offset eqn gives x=0.5.
        Site y0 = new LineSite(new Point(0, 0), new Point(1, 0), +1);
        Site y1 = new LineSite(new Point(0, 1), new Point(1, 1), +1);
        Site x0 = new LineSite(new Point(0, 0), new Point(0, 1), -1); // x=0, use k3=-1 for positive x solution

        LLLPARASolver solver = new LLLPARASolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(y0, y0.k(), y1, y1.k(), x0, x0.k(), out);

        assertEquals(1, n);
        Solution s = out.get(0);
        assertEquals(0.5, s.p.x, 1e-12);
        assertEquals(0.5, s.p.y, 1e-12);
        assertEquals(0.5, s.t, 1e-12);
    }
}

