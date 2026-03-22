package org.rogach.jopenvoronoi.geometry;

import static org.rogach.jopenvoronoi.util.Numeric.chop;
import static org.rogach.jopenvoronoi.util.Numeric.sq;

import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;

import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.util.Numeric;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Half-edge in the Voronoi diagram.
 * <p>
 * Besides the topological links ({@link #twin}, {@link #next}, {@link #prev}),
 * an edge stores an eight-parameter geometric parametrization for its x/y
 * coordinates:
 *
 * <pre>
 * x = x1 - x2 - x3*t +/- x4 * sqrt(square(x5 + x6*t) - square(x7 + x8*t))
 * </pre>
 *
 * <p>
 * The same formula is used for the y-coordinate. Depending on the adjacent
 * sites, this parametrization describes linear bisectors (line/line), parabolas
 * (point/line), or conics such as hyperbolas/ellipses (point/arc).
 */
public class Edge {

	public Vertex source;
	public Vertex target;
	public Edge twin;
	public Edge base;
	public Edge next;
	public Edge prev;
	public Face face;
	public Face nullFace;
	public boolean hasNullFace;
	/** offset-direction from the adjacent site, either +1 or -1 */
	public double k;
	public EdgeType type;
	public boolean valid;
	/** 8-parameter parametrization */
	private double[] x = new double[8];
	/** 8-parameter parametrization */
	private double[] y = new double[8];
	/** flag to choose either +/- in front of sqrt() */
	private boolean sign;
	/** true if ::LINESITE-edge inserted in this direction */
	public boolean insertedDirection;
	public int diagramIndex = -1;
	public int sourceOutIndex = -1;
	public int targetInIndex = -1;

	public Edge(Vertex source, Vertex target) {
		hasNullFace = false;
		valid = true;

		this.source = source;
		this.target = target;
		this.base = this;
	}

	public void copyFrom(Edge other) {
		this.sign = other.sign;
		this.face = other.face;
		this.k = other.k;
		this.nullFace = other.nullFace;
		this.hasNullFace = other.hasNullFace;
		this.type = other.type;
		this.valid = other.valid;
		this.insertedDirection = other.insertedDirection;
		x[0] = other.x[0];
		x[1] = other.x[1];
		x[2] = other.x[2];
		x[3] = other.x[3];
		x[4] = other.x[4];
		x[5] = other.x[5];
		x[6] = other.x[6];
		x[7] = other.x[7];
		y[0] = other.y[0];
		y[1] = other.y[1];
		y[2] = other.y[2];
		y[3] = other.y[3];
		y[4] = other.y[4];
		y[5] = other.y[5];
		y[6] = other.y[6];
		y[7] = other.y[7];
	}

	/**
	 * Return the point on this edge at the given offset-distance {@code t}.
	 * <p>
	 * The eight-parameter formula for a point on the edge is:
	 * 
	 * <pre>{@code
	 * x = x1 - x2 - x3*t +/- x4 * sqrt(square(x5+x6*t) - square(x7+x8*t))
	 * }</pre>
	 *
	 * @param t offset distance
	 * @return point on the edge at offset distance {@code t}
	 */
	public Point point(double t) {
		var discr1 = chop(sq(x[4] + x[5] * t) - sq(x[6] + x[7] * t), Numeric.STRICT_ZERO_EPSILON);
		var discr2 = chop(sq(y[4] + y[5] * t) - sq(y[6] + y[7] * t), Numeric.STRICT_ZERO_EPSILON);
		if ((discr1 >= 0) && (discr2 >= 0)) {
			double psig = sign ? +1 : -1;
			double nsig = sign ? -1 : +1;
			var xc = x[0] - x[1] - x[2] * t + psig * x[3] * Math.sqrt(discr1);
			var yc = y[0] - y[1] - y[2] * t + nsig * y[3] * Math.sqrt(discr2);
			if (xc != xc) { // test for NaN!
				throw new RuntimeException();
			}
			return new Point(xc, yc);
		} else {
			return new Point(x[0] - x[1] - x[2] * t, y[0] - y[1] - y[2] * t); // coordinates without sqrt()
		}
	}

	/**
	 * Returns the midpoint of the edge source and target positions
	 * 
	 * @return the midpoint between the source and target vertices
	 */
	public Point position() {
		return (source.position.add(target.position)).mult(0.5); // TODO pre-calculate?
	}

	/**
	 * Initialize the geometric parametrization for the bisector between two sites.
	 *
	 * @param s1  first adjacent site
	 * @param s2  second adjacent site
	 * @param sig selects the sign used in front of the square-root term
	 */
	public void setParameters(Site s1, Site s2, boolean sig) {
		sign = sig; // sqrt() sign for edge-parametrization
		if (s1.isPoint() && s2.isPoint()) {
			setPpParameters(s1, s2);
		} else if (s1.isPoint() && s2.isLine()) {
			setPlParameters(s1, s2);
		} else if (s2.isPoint() && s1.isLine()) { // LP
			setPlParameters(s2, s1);
			sign = !sign;
		} else if (s1.isLine() && s2.isLine()) {
			setLlParameters(s2, s1);
		} else if (s1.isPoint() && s2.isArc()) {
			setPaParameters(s1, s2);
		} else if (s2.isPoint() && s1.isArc()) { // AP
			sign = !sign;
			setPaParameters(s2, s1);

		} else if (s1.isLine() && s2.isArc()) {
			setLaParameters(s1, s2);
		} else if (s2.isLine() && s1.isArc()) {
			setLaParameters(s2, s1);
		} else {
			throw new RuntimeException("Unexpected combination of sites");
			// AA
		}
	}

	// Set parameters for a point/point bisector.
	private void setPpParameters(Site s1, Site s2) {
		assert (s1.isPoint() && s2.isPoint()) : " s1.isPoint() && s2.isPoint() ";
		var d = (s1.position().sub(s2.position())).norm();
		var alfa1 = (s2.x() - s1.x()) / d;
		var alfa2 = (s2.y() - s1.y()) / d;
		var alfa3 = -d / 2;

		type = EdgeType.LINE;
		x[0] = s1.x();
		x[1] = alfa1 * alfa3; //
		x[2] = 0;
		x[3] = -alfa2;
		x[4] = 0;
		x[5] = +1;
		x[6] = alfa3;
		x[7] = 0;
		y[0] = s1.y();
		y[1] = alfa2 * alfa3;
		y[2] = 0;
		y[3] = -alfa1;
		y[4] = 0;
		y[5] = +1;
		y[6] = alfa3;
		y[7] = 0;
	}

	// Set parameters for a point/line parabola bisector.
	private void setPlParameters(Site s1, Site s2) {
		assert (s1.isPoint() && s2.isLine()) : " s1.isPoint() && s2.isLine() ";

		type = EdgeType.PARABOLA;
		var alfa3 = s2.a() * s1.x() + s2.b() * s1.y() + s2.c(); // signed distance to line

		// figure out kk, i.e. offset-direction for LineSite
		x[0] = s1.x(); // xc1
		x[1] = s2.a() * alfa3; // alfa1*alfa3
		x[2] = s2.a(); // *kk; // -alfa1 = - a2 * k2?
		x[3] = s2.b(); // alfa2 = b2
		x[4] = 0; // alfa4 = r1 (PointSite has zero radius)
		x[5] = +1; // lambda1 (allways positive offset from PointSite)
		x[6] = alfa3; // alfa3= a2*xc1+b2*yc1+d2?
		x[7] = +1; // kk; // -1 = k2 side of line??

		y[0] = s1.y(); // yc1
		y[1] = s2.b() * alfa3; // alfa2*alfa3
		y[2] = s2.b(); // *kk; // -alfa2 = -b2
		y[3] = s2.a(); // alfa1 = a2
		y[4] = 0; // alfa4 = r1 (PointSite has zero radius)
		y[5] = +1; // lambda1 (allways positive offset from PointSite)
		y[6] = alfa3; // alfa3
		y[7] = +1; // kk; // -1 = k2 side of line??
	}

	/** Set parameters for a separator edge. */
	public void setSepParameters(Point endp, Point p) {
		type = EdgeType.SEPARATOR;
		var dx = p.x - endp.x;
		var dy = p.y - endp.y;
		var d = p.sub(endp).norm();
		x[0] = endp.x;
		x[2] = -dx / d; // negative of normalized direction from endp to p
		y[0] = endp.y;
		y[2] = -dy / d;

		x[1] = 0;
		x[3] = 0;
		x[4] = 0;
		x[5] = 0;
		x[6] = 0;
		x[7] = 0;
		y[1] = 0;
		y[3] = 0;
		y[4] = 0;
		y[5] = 0;
		y[6] = 0;
		y[7] = 0;
	}

	// Set parameters for a parallel line/line bisector.
	private void setLlParaParameters(Site s1, Site s2) {
		assert (s1.isLine() && s2.isLine()) : " s1.isLine() && s2.isLine() ";
		type = EdgeType.PARA_LINELINE;

		// find a point (x1,y1) on the line s1
		// ax+by+c=0
		var x1 = 0D;
		var y1 = 0D;
		if (Math.abs(s1.a()) > Math.abs(s1.b())) {
			y1 = 0;
			x1 = -s1.c() / s1.a();
		} else {
			x1 = 0;
			y1 = -s1.c() / s1.b();
		}

		// find a point (x2,y2) on the line s2
		// ax+by+c=0
		var x2 = 0D;
		var y2 = 0D;
		if (Math.abs(s2.a()) > Math.abs(s2.b())) {
			y2 = 0;
			x2 = -s2.c() / s2.a();
		} else {
			x2 = 0;
			y2 = -s2.c() / s2.b();
		}

		// now e.g. the s2 line is given by
		// p = (x2,y2) + t*(-b2, a)
		// and we can find the projection of (x1,y1) onto s2 as
		// p1 = p2 = p0 + t*v
		var p1 = new Point(x1, y1);
		var p2 = new Point(x2, y2);
		var v = new Point(-s2.b(), s2.a());
		var t = p1.sub(p2).dot(v) / v.dot(v);
		var p1_proj = p2.add(v.mult(t));

		assert (p1.sub(p1_proj).norm() > 0) : " p1.sub(p1_proj).norm() > 0 ";

		// from this point, go a distance d/2 in the direction of the normal
		// to find a point through which the bisector passes
		x1 = x1 + p1_proj.sub(p1).x / 2;
		y1 = y1 + p1_proj.sub(p1).y / 2;
		// the tangent of the bisector (as well as the two line-sites) is a vector
		// (-b , a)

		x[0] = x1;
		x[1] = -s1.b();
		y[0] = y1;
		y[1] = s1.a();

		x[2] = 0;
		x[3] = 0;
		x[4] = 0;
		x[5] = 0;
		x[6] = 0;
		x[7] = 0;
		y[2] = 0;
		y[3] = 0;
		y[4] = 0;
		y[5] = 0;
		y[6] = 0;
		y[7] = 0;
	}

	// Set parameters for a non-parallel line/line bisector.
	private void setLlParameters(Site s1, Site s2) { // Held thesis p96
		assert (s1.isLine() && s2.isLine()) : " s1.isLine() && s2.isLine() ";
		type = EdgeType.LINELINE;
		var delta = s1.a() * s2.b() - s1.b() * s2.a();

		// (numerically) parallel line segments - the generic LLL solver
		// is numerically unstable for parallel cases
		if (Math.abs(delta) <= Numeric.STRICT_ZERO_EPSILON) {
			setLlParaParameters(s1, s2);
			return;
		}

		assert (delta != 0) : " delta != 0 ";
		var alfa1 = (s1.b() * s2.c() - s2.b() * s1.c()) / delta;
		var alfa2 = (s2.a() * s1.c() - s1.a() * s2.c()) / delta;
		var alfa3 = -(s2.b() - s1.b()) / delta;
		var alfa4 = -(s1.a() - s2.a()) / delta;

		// point (alfa1,alfa2) is the intersection point between the line-segments
		// vector (-alfa3,-alfa4) is the direction/tangent of the bisector
		x[0] = alfa1;
		x[2] = -alfa3;
		y[0] = alfa2;
		y[2] = -alfa4;

		x[1] = 0;
		x[3] = 0;
		x[4] = 0;
		x[5] = 0;
		x[6] = 0;
		x[7] = 0;
		y[1] = 0;
		y[3] = 0;
		y[4] = 0;
		y[5] = 0;
		y[6] = 0;
		y[7] = 0;
	}

	// Set parameters for a point/arc bisector.
	private void setPaParameters(Site s1, Site s2) {
		assert (s1.isPoint() && s2.isArc()) : " s1.isPoint() && s2.isArc() ";
		type = EdgeType.HYPERBOLA; // hyperbola or ellipse?
		var lamb2 = 1.0;

		// distance between centers
		var d = Math.sqrt((s1.x() - s2.x()) * (s1.x() - s2.x()) + (s1.y() - s2.y()) * (s1.y() - s2.y()));
		assert (d > 0) : " d > 0 ";
		if (d <= s2.r()) {
			lamb2 = -1.0;
			sign = !sign;
		}

		var alfa1 = (s2.x() - s1.x()) / d;
		var alfa2 = (s2.y() - s1.y()) / d;
		var alfa3 = (s2.r() * s2.r() - d * d) / (2 * d);
		var alfa4 = (lamb2 * s2.r()) / d;
		x[0] = s1.x();
		x[1] = alfa1 * alfa3;
		x[2] = alfa1 * alfa4;
		x[3] = alfa2;
		x[4] = 0; // r1; PointSite has zero radius
		x[5] = +1; // lamb1; allways outward offset from PointSite
		x[6] = alfa3;
		x[7] = alfa4;

		y[0] = s1.y();
		y[1] = alfa2 * alfa3;
		y[2] = alfa2 * alfa4;
		y[3] = alfa1;
		y[4] = 0; // r1; PointSite has zero radius
		y[5] = +1; // lamb1; allways outward offset from PointSite
		y[6] = alfa3;
		y[7] = alfa4;
	}

	// Set parameters for a line/arc bisector.
	private void setLaParameters(Site s1, Site s2) {
		assert (s1.isLine() && s2.isArc()) : " s1.isLine() && s2.isArc() ";
		type = EdgeType.PARABOLA;
		double lamb2;
		if (s2.cw()) {
			lamb2 = +1.0;
		} else {
			lamb2 = -1.0;
		}
		var alfa1 = s1.a(); // a2
		var alfa2 = s1.b(); // b2
		var alfa3 = (s1.a() * s2.x() + s1.b() * s2.y() + s1.c());
		var alfa4 = s2.r();
		double kk = +1;

		x[0] = s2.x();
		x[1] = alfa1 * alfa3;
		x[2] = alfa1 * kk;
		x[3] = alfa2;
		x[4] = alfa4;
		x[5] = lamb2;
		x[6] = alfa3;
		x[7] = kk;

		y[0] = s2.y();
		y[1] = alfa2 * alfa3;
		y[2] = alfa2 * kk;
		y[3] = alfa1;
		y[4] = alfa4;
		y[5] = lamb2;
		y[6] = alfa3;
		y[7] = kk;
	}

	/**
	 * Returns the minimum {@code t}-value for this edge.
	 *
	 * @param s1 first site adjacent to this edge
	 * @param s2 second site adjacent to this edge
	 * @return minimum {@code t}-value for this edge. This function dispatches to a
	 *         helper function based on the sites {@code s1} and {@code s2}. It is
	 *         used only for positioning {@code APEX} vertices.
	 */
	public double minimumT(Site s1, Site s2) {
		if (s1.isPoint() && s2.isPoint()) {
			return minimumPpT(s1, s2);
		} else if (s1.isPoint() && s2.isLine()) {
			return minimumPlT(s1, s2);
		} else if (s2.isPoint() && s1.isLine()) {
			return minimumPlT(s2, s1);
		} else if (s1.isLine() && s2.isLine()) {
			return 0;
		} else if (s1.isPoint() && s2.isArc()) {
			return minimumPaT(s1, s2);
		} else if (s2.isPoint() && s1.isArc()) {
			return minimumPaT(s2, s1);
		} else {
			throw new RuntimeException("Unexpected site types");
			// todo: AP, AL, AA
		}
	}

	// Minimum t-value for a point/point edge.
	private double minimumPpT(Site s1, Site s2) {
		assert (s1.isPoint() && s2.isPoint()) : " s1.isPoint() && s2.isPoint() ";
		var p1p2 = s1.position().sub(s2.position()).norm();
		assert (p1p2 >= 0) : " p1p2 >=0 ";
		return p1p2 / 2; // this splits point-point edges at APEX
	}

	// Minimum t-value for a point/line parabola.
	private double minimumPlT(Site s1, Site s2) {
		var mint = -x[6] / (2.0 * x[7]);
		assert (mint >= 0) : " mint >=0 ";
		return mint;
	}

	// Minimum t-value for a point/arc edge.
	private double minimumPaT(Site s1, Site s2) {
		assert (s1.isPoint() && s2.isArc()) : " s1.isPoint() && s2.isArc() ";
		var p1p2 = s1.position().sub(s2.apexPoint(s1.position())).norm(); // - s2->r() ;
		assert (p1p2 >= 0) : " p1p2 >=0 ";
		return p1p2 / 2; // this splits point-point edges at APEX
	}

	/**
	 * Sample a MIC candidate point and its exact clearance radius at normalized
	 * parameter {@code u} along this edge.
	 *
	 * @param u normalized edge parameter from {@code 0.0} at {@link #source} to
	 *          {@code 1.0} at {@link #target}
	 * @return a pair containing the sampled point and its exact clearance radius
	 */
	public Entry<Point, Double> micSample(double u) {
		Point p;
		double r;

		if (type == EdgeType.LINE || type == EdgeType.LINELINE || type == EdgeType.PARA_LINELINE) {
			// Linear interpolation of the position.
			Point src = source.position;
			Point trg = target.position;
			p = src.add(trg.sub(src).mult(u));

			// Exact geometric radius: distance from p to the generating site.
			Site s = face.getSite();
			Point pa = s.apexPoint(p);
			r = p.sub(pa).norm();
		} else if (type == EdgeType.PARABOLA) {
			// For parabolas, the distance (radius) is the parameter t.
			// Interpolate t between the source and target clearance radii.
			double srcT = source.dist();
			double trgT = target.dist();
			double uT = srcT + u * (trgT - srcT);

			// Recover the position using the existing point(t) parametrization.
			p = point(uT);
			r = uT;
		} else {
			throw new RuntimeException("Unsupported edge type for micSample: " + type);
		}

		return new SimpleEntry<>(p, r);
	}

	/**
	 * Returns the two boundary touchpoints of the MIC sampled at normalized
	 * parameter {@code u} along this edge.
	 * <p>
	 * For a Voronoi / medial-axis edge, the sampled MIC center is equidistant from
	 * the two generator sites adjacent to the edge. The MIC touchpoints are the
	 * closest / support points on those two sites.
	 *
	 * @param u normalized edge parameter from {@code 0.0} at {@link #source} to
	 *          {@code 1.0} at {@link #target}
	 * @return array of length 2: {footA, footB}
	 */
	public Point[] micTouchPoints(double u) {
		Entry<Point, Double> sample = micSample(u);
		Point p = sample.getKey();

		Site s1 = adjacentSiteA();
		Site s2 = adjacentSiteB();

		if (s1 == null || s2 == null) {
			throw new RuntimeException("micFootPoints requires two adjacent sites on edge " + this);
		}

		Point f1 = s1.apexPoint(p);
		Point f2 = s2.apexPoint(p);

		return new Point[] { f1, f2 };
	}

	/**
	 * Returns {@code n} evenly-spaced sample points along this edge, including the
	 * start ({@link #source}) and end ({@link #target}) positions.
	 * <p>
	 * Sampling is performed by interpolating the offset-distance parameter
	 * {@code t} uniformly between {@code source.dist()} and {@code target.dist()},
	 * then evaluating {@link #point(double)} at each step.
	 *
	 * @param n number of sample points; must be at least 2
	 * @return list of {@code n} points sampled uniformly from source to target
	 * @throws IllegalArgumentException if {@code n} is less than 2
	 */
	public List<Point> samplePoints(int n) {
		if (n < 2) {
			throw new IllegalArgumentException("n must be at least 2, got " + n);
		}
		double srcT = source.dist();
		double trgT = target.dist();
		List<Point> points = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			double u = (double) i / (n - 1);
			double t = srcT + u * (trgT - srcT);
			points.add(point(t));
		}
		return points;
	}

	/**
	 * Computes the arc length of this edge analytically where possible. Returns the
	 * exact length for line segments and parabolas, and falls back to numerical
	 * approximation for higher-order conics like ellipses or hyperbolas.
	 * 
	 * @return the total arc length from the source to the target vertex
	 */
	public double length() {
		if (type == EdgeType.LINE || type == EdgeType.LINELINE || type == EdgeType.PARA_LINELINE || type == EdgeType.SEPARATOR || type == EdgeType.OUTEDGE
				|| type == EdgeType.LINESITE) {
			return source.position.distance(target.position);
		} else if (type == EdgeType.NULLEDGE) {
			return 0.0;
		} else if (type == EdgeType.PARABOLA) {
			// Extract focal length p algebraically using the stored 8-parameter edge
			// parametrization
			// $p$ natively equals $|x_4 x_5 - x_6 x_7| / 2$.
			double p = Math.abs(x[4] * x[5] - x[6] * x[7]) / 2.0;

			if (p > 1e-9) {
				// The parameter t is naturally the clearance radius (dist) at the vertex
				double t1 = source.dist();
				double t2 = target.dist();

				// Calculate the lateral offsets squared (X^2)
				double discr1 = sq(x[4] + x[5] * t1) - sq(x[6] + x[7] * t1);
				double discr2 = sq(x[4] + x[5] * t2) - sq(x[6] + x[7] * t2);

				// Use the lengths from the apex directly without branching
				double X1 = Math.sqrt(Math.max(0.0, discr1));
				double X2 = Math.sqrt(Math.max(0.0, discr2));

				return Math.abs(parabolaArcLength(X1, p) - parabolaArcLength(X2, p));
			} else {
				return source.position.distance(target.position);
			}
		}

		// Fallback for HYPERBOLA, ELLIPSE, or unhandled cases: Evaluate length via
		// adaptive recursive subdivision
		return lengthAdaptive(1e-5);
	}

	/**
	 * Fallback numerical arc length evaluating via adaptive recursive subdivision
	 * bounded by a deviation error, providing exceptionally accurate length where
	 * analytical forms are unavailable (like ellipses/hyperbolas).
	 */
	private double lengthAdaptive(double maxDev) {
		return lengthRecursive(source.dist(), target.dist(), source.position, target.position, maxDev * maxDev);
	}

	private double lengthRecursive(double t1, double t2, Point p1, Point p2, double maxDevSq) {
		double tMid = 0.5 * (t1 + t2);
		Point pMid = point(tMid);

		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		double chordLenSq = dx * dx + dy * dy;

		double deviationSq;
		if (chordLenSq < 1e-18) {
			double d = pMid.distance(p1);
			deviationSq = d * d;
		} else {
			double cross = dx * (p1.y - pMid.y) - (p1.x - pMid.x) * dy;
			deviationSq = (cross * cross) / chordLenSq;
		}

		if (deviationSq > maxDevSq && Math.abs(t2 - t1) > 1e-9) {
			return lengthRecursive(t1, tMid, p1, pMid, maxDevSq) + lengthRecursive(tMid, t2, pMid, p2, maxDevSq);
		} else {
			return Math.sqrt(chordLenSq);
		}
	}

	private double parabolaArcLength(double X, double p) {
		double u = X / (2 * p);
		double hyp = Math.sqrt(1 + u * u);
		double ln = (u >= 0) ? Math.log(u + hyp) : -Math.log(hyp - u);
		return p * (u * hyp + ln);
	}

	/**
	 * Returns adaptively sampled points along this edge, guaranteed to not deviate
	 * from the true piecewise curve by more than {@code maxDeviation}.
	 * 
	 * @param maxDeviation maximum allowed distance between the chord and true curve
	 * @return list of adaptively sampled points dynamically concentrated in highly
	 *         curved parts
	 */
	public List<Point> samplePoints(double maxDeviation) {
		if (maxDeviation <= 0) {
			throw new IllegalArgumentException("maxDeviation must be > 0");
		}

		if (type == EdgeType.LINE || type == EdgeType.LINELINE || type == EdgeType.PARA_LINELINE || type == EdgeType.SEPARATOR || type == EdgeType.OUTEDGE
				|| type == EdgeType.LINESITE || type == EdgeType.NULLEDGE) {
			List<Point> points = new ArrayList<>(2);
			points.add(source.position);
			if (type != EdgeType.NULLEDGE) {
				points.add(target.position);
			}
			return points;
		}

		List<Point> points = new ArrayList<>();
		points.add(source.position);
		samplePointsRecursive(source.dist(), target.dist(), source.position, target.position, maxDeviation * maxDeviation, points);
		points.add(target.position);
		return points;
	}

	private void samplePointsRecursive(double t1, double t2, Point p1, Point p2, double maxDevSq, List<Point> points) {
		double tMid = 0.5 * (t1 + t2);
		Point pMid = point(tMid);

		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		double chordLenSq = dx * dx + dy * dy;

		double deviationSq;
		if (chordLenSq < 1e-18) {
			double d = pMid.distance(p1);
			deviationSq = d * d;
		} else {
			double cross = dx * (p1.y - pMid.y) - (p1.x - pMid.x) * dy;
			deviationSq = (cross * cross) / chordLenSq;
		}

		if (deviationSq > maxDevSq && Math.abs(t2 - t1) > 1e-9) {
			samplePointsRecursive(t1, tMid, p1, pMid, maxDevSq, points);
			points.add(pMid);
			samplePointsRecursive(tMid, t2, pMid, p2, maxDevSq, points);
		}
	}

	/** First site adjacent to this edge. */
	private Site adjacentSiteA() {
		return (face != null) ? face.getSite() : null;
	}

	/** Second site adjacent to this edge. */
	private Site adjacentSiteB() {
		if (twin != null && twin.face != null) {
			return twin.face.getSite();
		}
		if (hasNullFace && nullFace != null) {
			return nullFace.getSite();
		}
		return null;
	}

	/**
	 * Returns true for non-null, non-null-face faces that have an attached site.
	 */
	private boolean isUsableFace(Face f) {
		return f != null && !f.isNullFace() && f.getSite() != null;
	}

	/** Returns true if the face is a usable point-site face. */
	private boolean isVertexFace(Face f) {
		return isUsableFace(f) && f.getSite().isPoint();
	}

	/** Returns true if the face is a usable line-segment-site face. */
	private boolean isSegmentFace(Face f) {
		return isUsableFace(f) && f.getSite().isLine();
	}

	/** Returns true if this edge separates two line-segment-site faces. */
	public boolean isLineLineBisector() {
		return twin != null && isSegmentFace(face) && isSegmentFace(twin.face);
	}

	/**
	 * Returns true if this edge separates a line-segment-site face and a point-site
	 * face.
	 */
	public boolean isLinePointBisector() {
		return twin != null && ((isSegmentFace(face) && isVertexFace(twin.face)) || (isVertexFace(face) && isSegmentFace(twin.face)));
	}

	/** Returns true if this edge separates two point-site faces. */
	public boolean isPointPointBisector() {
		return twin != null && isVertexFace(face) && isVertexFace(twin.face);
	}

	@Override
	public String toString() {
		return String.format("E(%s>%s)", source.position, target.position);
	}
}
