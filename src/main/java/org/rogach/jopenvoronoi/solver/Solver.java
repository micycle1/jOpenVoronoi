package org.rogach.jopenvoronoi.solver;

import java.util.List;

import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Solution;

/**
 * Abstract base class for Voronoi vertex position solvers.
 * 
 * The input to the solver is three {@link Site Sites} ({@code s1}, {@code s2},
 * {@code s3}) and three offset directions ({@code k1}, {@code k2},
 * {@code k3}). The output is a vector with one or more {@link Solution}
 * instances.
 */
public abstract class Solver {

	/**
	 * Solve for the position of a VoronoiVertex with the given adjacent sites and
	 * directions.
	 *
	 * @param s1 first adjacent Site
	 * @param k1 direction from {@code s1} to new VoronoiVertex
	 * @param s2 second adjacent Site
	 * @param k2 direction from {@code s2} to new VoronoiVertex
	 * @param s3 third adjacent Site
	 * @param k3 direction from {@code s3} to new VoronoiVertex
	 * @param slns Solution vector, will be updated by Solver
	 * @return number of solutions found
	 */
	public abstract int solve(Site s1, double k1, Site s2, double k2, Site s3, double k3, List<Solution> slns);

	// used by alt_sep_solver
	public void setType(int t) {
		type = t;
	}

	// flag for debug output
	// separator case type.
	// - type = 0 means l3 / p1 form a separator
	// - type = 1 means l3 / p2 form a separator
	int type;
}
