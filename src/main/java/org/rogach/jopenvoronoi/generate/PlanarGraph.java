package org.rogach.jopenvoronoi.generate;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class PlanarGraph {

	public List<Point2D> points;
	public List<Segment> segments;

	public PlanarGraph(List<Point2D> points, List<Segment> segments) {
		this.points = points;
		this.segments = segments;
	}

	public static class Segment {
		public Point2D stt;
		public Point2D end;

		public Segment(Point2D end, Point2D stt) {
			this.end = end;
			this.stt = stt;
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) {
				return true;
			}
			if (o == null || getClass() != o.getClass()) {
				return false;
			}

			var segment = (Segment) o;

			if (end != null ? !end.equals(segment.end) : segment.end != null) {
				return false;
			}
			if (stt != null ? !stt.equals(segment.stt) : segment.stt != null) {
				return false;
			}

			return true;
		}

		@Override
		public int hashCode() {
			var result = stt != null ? stt.hashCode() : 0;
			result = 31 * result + (end != null ? end.hashCode() : 0);
			return result;
		}
	}

	public static PlanarGraph fromPolygon(List<Point2D> points) {
		List<Segment> segments = new ArrayList<>();
		for (var q = 0; q < points.size(); q++) {
			if (q != points.size() - 1) {
				segments.add(new Segment(points.get(q), points.get(q + 1)));
			} else {
				segments.add(new Segment(points.get(points.size() - 1), points.get(0)));
			}
		}
		return new PlanarGraph(points, segments);
	}

	public VoronoiDiagram buildVoronoiDiagram() {
		Map<Point2D, Vertex> vertices = new HashMap<>();
		var vd = new VoronoiDiagram();
		for (Point2D p : points) {
			vertices.put(p, vd.insert_point_site(new Point(p.getX(), p.getY())));
		}
		for (Segment s : segments) {
			vd.insert_line_site(vertices.get(s.stt), vertices.get(s.end));
		}
		return vd;
	}

	public void buildIntoVoronoiDiagram(VoronoiDiagram v) {
		Map<Point2D, Vertex> vertices = new HashMap<>();
		for (Point2D p : points) {
			vertices.put(p, v.insert_point_site(new Point(p.getX(), p.getY())));
		}
		for (Segment s : segments) {
			v.insert_line_site(vertices.get(s.stt), vertices.get(s.end));
		}
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}

		var that = (PlanarGraph) o;

		if (points != null ? !points.equals(that.points) : that.points != null) {
			return false;
		}
		if (segments != null ? !segments.equals(that.segments) : that.segments != null) {
			return false;
		}

		return true;
	}

	@Override
	public int hashCode() {
		var result = points != null ? points.hashCode() : 0;
		result = 31 * result + (segments != null ? segments.hashCode() : 0);
		return result;
	}

	public static PlanarGraph minimizeFailure(PlanarGraph input) {
		System.out.printf("Minimizing input with %d points and %d segments\n", input.points.size(),
				input.segments.size());
		Throwable origException = null;
		try {
			input.buildVoronoiDiagram();
		} catch (Throwable e) {
			origException = e;
		}
		if (origException == null) {
			System.out.println("No exception, nothing to minimize");
			return input;
		}
		PlanarGraph current = input;
		int batch = current.points.size() / 2;
		while (true) {
			System.out.printf("\nAt the start of the iteration: %d points and %d segments\n", current.points.size(),
					current.segments.size());
			int c = 0;
			for (batch = batch > 0 ? batch : 1; batch >= 1; batch /= 2) {
				c = 0;
				System.out.printf("@%dx%d@", batch, current.points.size() / batch);
				System.out.flush();
				for (int offset = 0; offset + batch <= current.segments.size(); offset += batch) {
					List<PlanarGraph.Segment> lessSegments = new ArrayList<>(current.segments);
					for (int q = offset; q < offset + batch; q++) {
						lessSegments.remove(current.segments.get(q));
					}
					PlanarGraph modified = new PlanarGraph(current.points, lessSegments);
					try {
						modified.buildVoronoiDiagram();
						System.out.printf("|");
						System.out.flush();
					} catch (Throwable t) {
						current = modified;
						for (int q = 0; q < batch; q++) {
							System.out.printf("-");
						}
						System.out.flush();
						c += batch;
					}
				}

				for (int offset = 0; offset + batch <= current.points.size(); offset += batch) {
					List<Point2D> lessPoints = new ArrayList<>(current.points);
					int pointsRemoved = 0;
					for (int q = offset; q < offset + batch; q++) {
						Point2D p = current.points.get(q);
						boolean includedInSegment = false;
						for (PlanarGraph.Segment s : current.segments) {
							if (s.stt.equals(p) || s.end.equals(p)) {
								includedInSegment = true;
								break;
							}
						}
						if (!includedInSegment) {
							lessPoints.remove(current.points.get(q));
							pointsRemoved++;
						}
					}
					if (pointsRemoved > 0) {
						PlanarGraph modified = new PlanarGraph(lessPoints, current.segments);
						try {
							modified.buildVoronoiDiagram();
							System.out.printf("*");
							System.out.flush();
						} catch (Throwable t) {
							current = modified;
							for (int q = 0; q < pointsRemoved; q++) {
								System.out.printf(".");
							}
							System.out.flush();
							c += pointsRemoved;
						}
					} else {
						System.out.printf("*");
					}
				}
			}
			if (c == 0) {
				break;
			}
		}
		System.out.printf("\nMinimized input to %d points and %d segments\n", current.points.size(),
				current.segments.size());
		return current;
	}

	public static boolean isSelfIntersected(PlanarGraph input) {
		// check for intersecting segments
		for (PlanarGraph.Segment s1 : input.segments) {
			for (PlanarGraph.Segment s2 : input.segments) {
				if (!s1.equals(s2) && // do not compare segment with itself
				// do not count connected segments as intersecting
						!s1.stt.equals(s2.end) && !s2.stt.equals(s1.end) && !s1.end.equals(s2.end)) {
					if (Line2D.linesIntersect(s1.stt.getX(), s1.stt.getY(), s1.end.getX(), s1.end.getY(), s2.stt.getX(),
							s2.stt.getY(), s2.end.getX(), s2.end.getY())) {
						return true;
					}
				}
			}
		}

		// check for points lying directly on other segments
		for (PlanarGraph.Segment s : input.segments) {
			for (Point2D p : input.points) {
				if (!p.equals(s.stt) && !p.equals(s.end)) {
					if (Line2D.ptSegDist(s.stt.getX(), s.stt.getY(), s.end.getX(), s.end.getY(), p.getX(),
							p.getY()) < 1e-10) {
						return true;
					}
				}
			}
		}

		return false;
	}
}