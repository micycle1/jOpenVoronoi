package org.rogach.jopenvoronoi.pocket;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * Maximal Inscribed Circle (MIC) — describes one machining step in medial-axis
 * pocketing.
 * <p>
 * Each MIC pairs a previously-cut circle ({@link #c1}/{@link #r1}) with a new
 * circle to cut ({@link #c2}/{@link #r2}), together with four bi-tangent points
 * ({@link #t1}–{@link #t4}) that define the material boundaries between the two
 * circles. A downstream algorithm (see {@link PocketPath}) converts the
 * sequence of MICs into a continuous toolpath of line and arc segments.
 * <p>
 * When {@link #newBranch} is {@code true}, the tool must retract (or rapid)
 * from the end of the previous branch to the branch origin stored in
 * {@link #cPrev}/{@link #rPrev} before cutting this MIC. This occurs at
 * junctions of the medial axis — for example, the center of a square pocket
 * where the medial axis forks into four branches.
 */
public class MIC {

	/** center of the previously cut circle */
	public Point c1;
	/** center of the new circle */
	public Point c2;
	/** radius of the previously cut circle */
	public double r1;
	/** radius of the new circle */
	public double r2;
	/** first bi-tangent point */
	public Point t1;
	/** second bi-tangent point */
	public Point t2;
	/** third bi-tangent point */
	public Point t3;
	/** fourth bi-tangent point */
	public Point t4;
	/** touch points 1 and 2 of the circle */
	public Point tp1a, tp1b;
	/** touch points 1 and 2 of the previous circle */
	public Point tp2a, tp2b;
	/** is this a new branch? */
	public boolean newBranch;
	/** for a new branch, the previous center */
	public Point cPrev;
	/** for a new branch, the previous radius */
	public double rPrev;
}
