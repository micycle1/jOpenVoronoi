package org.rogach.jopenvoronoi.site;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.offset.LineOfs;
import org.rogach.jopenvoronoi.offset.Ofs;

//line segment Site
public class LineSite extends Site {

	Point _start; // < start Point of LineSite
	Point _end; // < end Point of LineSite
	private Point _p; // < position
	public Edge e; // < edge_descriptor to the ::LINESITE pseudo-edge

	public LineSite(Point st, Point en, double koff) {
		this(st, en, koff, null);
	}

	public LineSite(Site s) {
		this.eq = s.eqp();
		this.face = s.face;
		this._start = s.start();
		this._end = s.end();
		this._p = s.position(); // inherit position
	}

	// create line-site between start and end Point.
	public LineSite(Point st, Point en, double koff, Face f) {
		this._start = st;
		this._end = en;
		this._p = new Point((start().add(end()).mult(0.5)));
		face = f;
		eq.q = false;
		eq.a = _end.y - _start.y;
		eq.b = _start.x - _end.x;
		eq.k = koff; // ??
		eq.c = _end.x * _start.y - _start.x * _end.y;
		// now normalize
		var d = Math.sqrt(eq.a * eq.a + eq.b * eq.b);
		eq.a /= d;
		eq.b /= d;
		eq.c /= d;
		assert (Math.abs(eq.a * eq.a + eq.b * eq.b - 1.0) < 1e-5) : " Math.abs( eq.a*eq.a + eq.b*eq.b -1.0 ) < 1e-5";
	}

	@Override
	public Ofs offset(Point p1, Point p2) {
		return new LineOfs(p1, p2);
	}

	// closest point on start-end segment to given point.
	// project onto line and return either the projected point
	// or one endpoint of the linesegment
	@Override
	public Point apex_point(Point p) {
		var s_p = p.sub(_start);
		var s_e = _end.sub(_start);
		var t = s_p.dot(s_e) / s_e.dot(s_e);
		if (t < 0) {
			return _start;
		}
		if (t > 1) {
			return _end;
		} else {
			return _start.add(_end.sub(_start).mult(t));
		}
	}

	@Override
	public boolean in_region(Point p) {
		var t = in_region_t(p);
		return ((t >= 0) && (t <= 1));
	}

	@Override
	public double in_region_t(Point p) {
		var s_p = p.sub(_start);
		var s_e = _end.sub(_start);
		var t = s_p.dot(s_e) / s_e.dot(s_e);
		var eps = 1e-7;
		if (Math.abs(t) < eps) {
			t = 0.0;
		} else if (Math.abs(t - 1.0) < eps) {
			t = 1.0;
		}
		return t;
	}

	@Override
	public double in_region_t_raw(Point p) {
		var s_p = p.sub(_start);
		var s_e = _end.sub(_start);
		var t = s_p.dot(s_e) / s_e.dot(s_e);
		return t;
	}

	@Override
	public boolean isLine() {
		return true;
	}

	@Override
	public Point position() {
		return _p;
	}

	@Override
	public double a() {
		return eq.a;
	}

	@Override
	public double b() {
		return eq.b;
	}

	@Override
	public double c() {
		return eq.c;
	}

	@Override
	public double k() {
		assert (eq.k == 1 || eq.k == -1) : " eq.k==1 || eq.k==-1 ";
		return eq.k;
	}

	@Override
	public void set_c(Point p) {
		eq.c = -(eq.a * p.x + eq.b * p.y);
	}

	@Override
	public Point start() {
		return _start;
	}

	@Override
	public Point end() {
		return _end;
	}

	@Override
	public Edge edge() {
		return e;
	}

	@Override
	public String toString() {
		return String.format("LS(%s>%s)", _start, _end);
	}
};
