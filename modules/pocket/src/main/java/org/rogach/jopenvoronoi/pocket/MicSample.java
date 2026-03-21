package org.rogach.jopenvoronoi.pocket;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * A sampled medial-axis circle: center/radius plus the two boundary footpoints
 * where the clearance circle touches the pocket boundary.
 */
public final class MicSample {
	public final Point center;
	public final double radius;
	public final Point footA;
	public final Point footB;

	public MicSample(Point center, double radius, Point footA, Point footB) {
		this.center = center;
		this.radius = radius;
		this.footA = footA;
		this.footB = footB;
	}

	public boolean hasFootpoints() {
		return footA != null && footB != null;
	}
}