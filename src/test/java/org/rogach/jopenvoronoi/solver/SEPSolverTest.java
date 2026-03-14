package org.rogach.jopenvoronoi.solver;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Tests for SEPSolver (separator edge: one line-site and its endpoint point-site).
 */
public class SEPSolverTest {

    @Test
    public void separatorWithThirdLine() {
        // s1: line y=0 (k irrelevant here), s2: point at origin (0,0), s3: line x=1 with k3=+1.
        // Separator direction for s1 (y=0) is (0, 1), so solution should be at (0,1) with t=1.
        Site s1 = new LineSite(new Point(0, 0), new Point(1, 0), +1);
        Site s2 = new PointSite(new Point(0, 0));
        Site s3 = new LineSite(new Point(1, -1), new Point(1, 1), +1); // x = 1

        SEPSolver solver = new SEPSolver();
        List<Solution> out = new ArrayList<>();
        int n = solver.solve(s1, s1.k(), s2, 1, s3, s3.k(), out);

        assertEquals(1, n);
        Solution sol = out.get(0);
        assertEquals(0.0, sol.p.x, 1e-12);
        assertEquals(1.0, sol.p.y, 1e-12);
        assertEquals(1.0, sol.t, 1e-12);
    }
}

