package org.rogach.jopenvoronoi.generate;

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
}