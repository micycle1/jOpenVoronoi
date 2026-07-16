package org.rogach.jopenvoronoi.solver;

import static org.rogach.jopenvoronoi.util.Numeric.chop;
import static org.rogach.jopenvoronoi.util.Numeric.determinant;

import java.util.ArrayList;
import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Line-line-line solver.
 * <p>
 * Solves a 3x3 system.
 */
public class LLLSolver extends Solver {

//  a1 x + b1 y + c1 + k1 t = 0
//  a2 x + b2 y + c2 + k2 t = 0
//  a3 x + b3 y + c3 + k3 t = 0
//
// or in matrix form
//
//  ( a1 b1 k1 ) ( x )    ( c1 )
//  ( a2 b2 k2 ) ( y ) = -( c2 )          Ax = b
//  ( a3 b3 k3 ) ( t )    ( c3 )
//
//  Cramers rule x_i = det(A_i)/det(A)
//  where A_i is A with column i replaced by b

	// reusable per-solve scratch; solver instances are used single-threaded
	// (see VertexPositioner.solutionsBuffer for the same reasoning)
	private final Eq[] eq = { new Eq(), new Eq(), new Eq() };
	private final Site[] sites = new Site[3];
	private final double[] kvals = new double[3];
	private final LLLPARASolver paraSolver = new LLLPARASolver();
	private final List<Solution> paraSolutions = new ArrayList<>();

	@Override
	public int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns) {

		assert (s1.isLine() && s2.isLine() && s3.isLine()) : " s1.isLine() && s2.isLine() && s3.isLine() ";

		sites[0] = s1;
		sites[1] = s2;
		sites[2] = s3;
		kvals[0] = k1;
		kvals[1] = k2;
		kvals[2] = k3;
		for (var i = 0; i < 3; i++) {
			sites[i].eqp(kvals[i], eq[i]); // equation-parameters, in quad-precision
		}

		int i = 0, j = 1, k = 2;
		var d = chop(determinant(eq[i].a, eq[i].b, eq[i].k, eq[j].a, eq[j].b, eq[j].k,
				eq[k].a, eq[k].b, eq[k].k));
		var det_eps = 1e-6;
		if (Math.abs(d) > det_eps) {
			var t = determinant(eq[i].a, eq[i].b, -eq[i].c, eq[j].a, eq[j].b, -eq[j].c,
					eq[k].a, eq[k].b, -eq[k].c) / d;
			if (t >= 0) {
				var sol_x = determinant(-eq[i].c, eq[i].b, eq[i].k, -eq[j].c, eq[j].b, eq[j].k,
						-eq[k].c, eq[k].b, eq[k].k) / d;
				var sol_y = determinant(eq[i].a, -eq[i].c, eq[i].k, eq[j].a, -eq[j].c, eq[j].k,
						eq[k].a, -eq[k].c, eq[k].k) / d;

				slns.add(new Solution(new Point(sol_x, sol_y), t, k3)); // kk3 just passes through without any effect!?
				return 1;
			}
		} else {
			// Try parallel solver as fallback, if the small determinant is due to nearly
			// parallel edges
			for (i = 0; i < 3; i++) {
				j = (i + 1) % 3;
				var delta = Math.abs(eq[i].a * eq[j].b - eq[j].a * eq[i].b);
				if (delta <= 1e-14) {
					paraSolutions.clear();
					paraSolver.solve(sites[i], kvals[i], sites[j], kvals[j], sites[(i + 2) % 3], kvals[(i + 2) % 3],
							paraSolutions);
					var solution_count = 0;
					for (Solution s : paraSolutions) {
						// check that solution has proper offset-direction
						if (s3.end().sub(s3.start()).cross(s.p.sub(s3.start())) * k3 >= 0) {
							slns.add(s);
							solution_count++;
						}
					}
					return solution_count;
				}
			}
		}
		return 0; // no solution if determinant zero, or t-value negative
	}

}
