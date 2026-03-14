package org.rogach.jopenvoronoi.site;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.offset.LineOfs;
import org.rogach.jopenvoronoi.offset.Ofs;
import org.rogach.jopenvoronoi.solver.Eq;

/**
 * Representation of a line-segment generator (LineSite) for the Voronoi
 * diagram.
 *
 * <p>
 * A {@code LineSite} models a finite line segment with a normalized line
 * equation {@code a*x + b*y + c = 0} (where {@code a*a + b*b = 1}) and an
 * offset-direction {@code k} which must be {@code +1} or {@code -1}. The offset
 * equation used by solvers is:
 *
 * <pre>
 * a*x + b*y + c + k * t = 0
 * </pre>
 *
 * <p>
 * Note: callers are expected to pass {@code koff} equal to {@code +1} or
 * {@code -1}; many algorithms rely on that invariant.
 */
public class LineSite extends Site {

	/** start Point of LineSite */
	Point _start;
	/** end Point of LineSite */
	Point _end;
	/**
	 * midpoint used as the representative position
	 */
	private Point _p;
	/**
	 * the pseudo-edge associated with this LineSite (if any)
	 */
	public Edge e;

	/**
	 * Create a LineSite from segment endpoints and an offset-direction.
	 *
	 * @param st   segment start point
	 * @param en   segment end point
	 * @param koff offset-direction, must be {@code +1} or {@code -1}
	 */
	public LineSite(Point st, Point en, double koff) {
		this(st, en, koff, null);
	}

	/**
	 * Copy-constructor: build a LineSite from an existing Site. The line equation
	 * parameters are copied from the source site.
	 *
	 * @param s site with compatible line parameters
	 */
	public LineSite(Site s) {
		this.eq = new Eq(s.eqp());
		this.face = s.face;
		this._start = s.start();
		this._end = s.end();
		this._p = s.position(); // inherit position
	}

	/**
	 * Create a LineSite from endpoints, offset-direction and optional face.
	 *
	 * <p>
	 * Constructor computes and normalizes the line equation parameters {@code a},
	 * {@code b} and {@code c} and stores {@code koff} into {@code eq.k}.
	 *
	 * @param st   start point
	 * @param en   end point
	 * @param koff offset-direction (+1 or -1)
	 * @param f    optional face associated with this site (may be {@code null})
	 */
	public LineSite(Point st, Point en, double koff, Face f) {
		this._start = st;
		this._end = en;
		this._p = new Point((start().add(end()).mult(0.5))); // midpoint
		face = f;
		eq.q = false;
		eq.a = _end.y - _start.y;
		eq.b = _start.x - _end.x;
		eq.k = koff; // offset direction: +1 or -1
		eq.c = _end.x * _start.y - _start.x * _end.y;
		// now normalize
		var d = Math.sqrt(eq.a * eq.a + eq.b * eq.b);
		eq.a /= d;
		eq.b /= d;
		eq.c /= d;
		assert (Math.abs(eq.a * eq.a + eq.b * eq.b - 1.0) < 1e-5) : " Math.abs( eq.a*eq.a + eq.b*eq.b -1.0 ) < 1e-5";
	}

	/**
	 * Return an {@link Ofs} implementation describing the offset between two points
	 * on this site (used by solver machinery).
	 */
	@Override
	public Ofs offset(Point p1, Point p2) {
		return new LineOfs(p1, p2);
	}

	/**
	 * Return the closest point on the finite segment to {@code p}. This projects
	 * onto the infinite line and clamps to segment endpoints.
	 *
	 * @param p query point
	 * @return closest point on the segment
	 */
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

	/**
	 * Test whether point {@code p} lies within the parametric region of the segment
	 * (inclusive).
	 *
	 * @param p query point
	 * @return {@code true} if the projection parameter t is in [0,1]
	 */
	@Override
	public boolean in_region(Point p) {
		var t = in_region_t(p);
		return ((t >= 0) && (t <= 1));
	}

	/**
	 * Return the parametric projection {@code t} (clamped) of {@code p} onto the
	 * segment, with numerical tolerance applied to near-0/1 values.
	 *
	 * @param p query point
	 * @return {@code t} in [0,1] (approx.)
	 */
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

	/**
	 * Offset direction parameter.
	 *
	 * @return {@code +1} or {@code -1} (asserted)
	 */
	@Override
	public double k() {
		assert (eq.k == 1 || eq.k == -1) : " eq.k==1 || eq.k==-1 ";
		return eq.k;
	}

	/**
	 * Set {@code c} such that the normalized line equation passes through point
	 * {@code p}.
	 *
	 * @param p point to set into the line equation
	 */
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
}
