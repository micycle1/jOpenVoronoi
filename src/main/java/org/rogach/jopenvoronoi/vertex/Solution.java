package org.rogach.jopenvoronoi.vertex;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * A new vertex position solution (position, offset-distance, side) includes the
 * offset-distance {@code t} and the offset direction {@code k3}.
 */
public class Solution {
	// position
	public Point p;
	// clearance-disk radius
	public double t;
	// offset direction to third adjacent Site
	public double k3;

	/**
	 * @param pt vertex position
	 * @param tv clearance-disk radius
	 * @param k offset direction
	 */
	public Solution(Point pt, double tv, double k) {
		this.p = pt;
		this.t = tv;
		this.k3 = k;
	}

	@Override
	public String toString() {
		return String.format("Solution(p = %s, t = %s, k = %d)", p, t, (int) k3);
	}
}
