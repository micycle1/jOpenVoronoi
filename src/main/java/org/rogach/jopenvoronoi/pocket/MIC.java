package org.rogach.jopenvoronoi.pocket;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * Maximal Inscribed Circle — a pair of clearance disks plus bi-tangent data.
 * <p>
 * It is the responsibility of a downstream algorithm to lay down a toolpath
 * that machines the area between {@link #c2}/{@link #r2} (the new circle) and
 * {@link #c1}/{@link #r1} (the circle that is assumed already cut).
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
	/** is this a new branch? */
	public boolean newBranch;
	/** for a new branch, the previous center */
	public Point cPrev;
	/** for a new branch, the previous radius */
	public double rPrev;
}
