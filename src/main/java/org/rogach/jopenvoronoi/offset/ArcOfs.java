package org.rogach.jopenvoronoi.offset;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * offset-element of PointSite or ArcSite
 */
public class ArcOfs extends Ofs {
	/** start */
	private Point _start;
	/** end */
	private Point _end;
	/** center */
	private Point c;
	/** radius */
	private double r;

	/**
	 * @param p1 start Point
	 * @param p2 end Point
	 * @param cen center Point
	 * @param rad radius
	 */
	public ArcOfs(Point p1, Point p2, Point cen, double rad) {
		this._start = p1;
		this._end = p2;
		this.c = cen;
		this.r = rad;
	}

	@Override
	public String toString() {
		return String.format("ArcOfs from %s to %s r=%f\n", _start, _end, r);
	}

	@Override
	public double radius() {
		return r;
	}

	@Override
	public Point center() {
		return c;
	}

	@Override
	public Point start() {
		return _start;
	}

	@Override
	public Point end() {
		return _end;
	}
};
