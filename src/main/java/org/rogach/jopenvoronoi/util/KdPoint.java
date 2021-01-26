package org.rogach.jopenvoronoi.util;

import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;

public class KdPoint {

	public Point p;
	public Face face;

	public KdPoint(Point p, Face face) {
		this.p = p;
		this.face = face;
	}
}
