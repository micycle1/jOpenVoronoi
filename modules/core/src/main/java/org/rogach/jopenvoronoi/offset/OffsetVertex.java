package org.rogach.jopenvoronoi.offset;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;

/**
 * Line- or arc-vertex of an offset curve. TODO this duplicates the idea of the Ofs class. Remove this or Ofs!
 */
public class OffsetVertex {
	/** position (start) */
	public Point p;
	/** arc radius (line-vertex is indicated by radius of -1) */
	public double r;
	/** arc center */
	public Point c;
	/** clockwise (or not) */
	public boolean cw;
	/** corresponding face in the vd-graph */
	public Face f;
	public Edge e; // corresponding edge

	public OffsetVertex(Point p, double r, Point c, boolean cw, Face f, Edge e) {
		this.p = p;
		this.r = r;
		this.c = c;
		this.cw = cw;
		this.f = f;
		this.e = e;
	}

	public OffsetVertex(Point p, Edge e) {
		this.p = p;
		r = -1;
		c = null;
		cw = false;
		f = null;
		this.e = e;
	}
}
