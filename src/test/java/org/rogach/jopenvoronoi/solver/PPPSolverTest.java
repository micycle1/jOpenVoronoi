package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for PPPSolver (three point-sites).
 */
public class PPPSolverTest {

    @Test
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

        double eps = 1e-9;
        assertEquals(0.5, s.p.x, 1e-9);
        assertEquals(0.5, s.p.y, 1e-9);
        assertEquals(Math.sqrt(0.5), s.t, 1e-9);
    }
}

