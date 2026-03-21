package org.rogach.jopenvoronoi.jts;

import java.util.ArrayList;
import java.util.Collection;
import java.util.IdentityHashMap;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Bridges JTS geometries and jOpenVoronoi's native geometry types.
 */
public final class JtsGeometryIO {

	private JtsGeometryIO() {
	}

	public static VoronoiDiagram toVoronoiDiagram(Geometry geometry) {
		if (geometry == null) {
			throw new IllegalArgumentException("Geometry cannot be null.");
		}

		VoronoiDiagram diagram = new VoronoiDiagram();
		addGeometry(diagram, geometry);
		return diagram;
	}

	static void addGeometry(VoronoiDiagram diagram, Geometry geometry) {
		if (diagram == null) {
			throw new IllegalArgumentException("Voronoi diagram cannot be null.");
		}
		if (geometry == null) {
			throw new IllegalArgumentException("Geometry cannot be null.");
		}

		List<Vertex> vertices = new ArrayList<>();
		addPointalGeometry(diagram, geometry, vertices);
		addLinearGeometry(diagram, geometry, vertices);
	}

	public static List<Vertex> addLineStringSites(VoronoiDiagram diagram, LineString lineString) {
		if (diagram == null) {
			throw new IllegalArgumentException("Voronoi diagram cannot be null.");
		}
		if (lineString == null) {
			throw new IllegalArgumentException("LineString cannot be null.");
		}

		Coordinate[] coordinates = lineString.getCoordinates();
		if (coordinates.length < 2) {
			throw new IllegalArgumentException("LineString must contain at least two coordinates.");
		}

		boolean closed = lineString.isClosed();
		int limit = closed && coordinates.length > 1 ? coordinates.length - 1 : coordinates.length;
		if (limit < 2) {
			throw new IllegalArgumentException("LineString must contain at least two distinct coordinates.");
		}

		List<Vertex> vertices = new ArrayList<>(limit);
		for (int i = 0; i < limit; i++) {
			Coordinate coordinate = coordinates[i];
			vertices.add(diagram.insertPointSite(toPoint(coordinate)));
		}

		for (int i = 0; i < vertices.size() - 1; i++) {
			diagram.insertLineSite(vertices.get(i), vertices.get(i + 1));
		}
		if (closed) {
			diagram.insertLineSite(vertices.get(vertices.size() - 1), vertices.get(0));
		}

		return vertices;
	}

	public static List<Vertex> addPolygonSites(VoronoiDiagram diagram, Polygon polygon) {
		if (polygon == null) {
			throw new IllegalArgumentException("Polygon cannot be null.");
		}

		List<Vertex> vertices = new ArrayList<>();
		vertices.addAll(addLinearRingSites(diagram, polygon.getExteriorRing()));
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			vertices.addAll(addLinearRingSites(diagram, polygon.getInteriorRingN(i)));
		}
		return vertices;
	}

	public static List<Vertex> addLinearRingSites(VoronoiDiagram diagram, LinearRing ring) {
		if (ring == null) {
			throw new IllegalArgumentException("LinearRing cannot be null.");
		}
		return addLineStringSites(diagram, ring);
	}

	public static MultiLineString toMultiLineString(VoronoiDiagram diagram) {
		return toMultiLineString(diagram, new GeometryFactory());
	}

	public static MultiLineString toMultiLineString(VoronoiDiagram diagram, GeometryFactory geometryFactory) {
		if (diagram == null) {
			throw new IllegalArgumentException("Voronoi diagram cannot be null.");
		}
		return toMultiLineString(diagram.getDiagram(), geometryFactory);
	}

	public static MultiLineString toMultiLineString(HalfEdgeDiagram diagram) {
		return toMultiLineString(diagram, new GeometryFactory());
	}

	public static MultiLineString toMultiLineString(HalfEdgeDiagram diagram, GeometryFactory geometryFactory) {
		if (diagram == null) {
			throw new IllegalArgumentException("HalfEdgeDiagram cannot be null.");
		}
		return toMultiLineString(diagram.edges, geometryFactory);
	}

	public static MultiLineString toMultiLineString(Collection<Edge> edges, GeometryFactory geometryFactory) {
		if (edges == null) {
			throw new IllegalArgumentException("Edge collection cannot be null.");
		}
		if (geometryFactory == null) {
			throw new IllegalArgumentException("GeometryFactory cannot be null.");
		}

		IdentityHashMap<Edge, Boolean> seenBases = new IdentityHashMap<>();
		List<LineString> lineStrings = new ArrayList<>();
		for (Edge edge : edges) {
			if (edge == null || edge.base == null || !edge.valid || edge.source == null || edge.target == null
					|| edge.source.position == null || edge.target.position == null) {
				continue;
			}
			if (seenBases.put(edge.base, Boolean.TRUE) != null) {
				continue;
			}
			lineStrings.add(geometryFactory.createLineString(
					new Coordinate[] { toCoordinate(edge.source.position), toCoordinate(edge.target.position) }));
		}
		return geometryFactory.createMultiLineString(lineStrings.toArray(new LineString[0]));
	}

	public static Coordinate toCoordinate(org.rogach.jopenvoronoi.geometry.Point point) {
		if (point == null) {
			throw new IllegalArgumentException("Point cannot be null.");
		}
		return new Coordinate(point.x, point.y);
	}

	public static org.rogach.jopenvoronoi.geometry.Point toPoint(Coordinate coordinate) {
		if (coordinate == null) {
			throw new IllegalArgumentException("Coordinate cannot be null.");
		}
		return new org.rogach.jopenvoronoi.geometry.Point(coordinate.getX(), coordinate.getY());
	}

	private static void addPointalGeometry(VoronoiDiagram diagram, Geometry geometry, List<Vertex> vertices) {
		if (geometry instanceof Point) {
			addPointSite(diagram, (Point) geometry, vertices);
			return;
		}
		if (geometry instanceof GeometryCollection) {
			GeometryCollection collection = (GeometryCollection) geometry;
			for (int i = 0; i < collection.getNumGeometries(); i++) {
				addPointalGeometry(diagram, collection.getGeometryN(i), vertices);
			}
			return;
		}
	}

	private static void addLinearGeometry(VoronoiDiagram diagram, Geometry geometry, List<Vertex> vertices) {
		if (geometry instanceof Point) {
			return;
		}
		if (geometry instanceof LineString) {
			vertices.addAll(addLineStringSites(diagram, normalizeClosedLineal((LineString) geometry)));
			return;
		}
		if (geometry instanceof Polygon) {
			vertices.addAll(addPolygonSites(diagram, normalizePolygon((Polygon) geometry)));
			return;
		}
		if (geometry instanceof GeometryCollection) {
			GeometryCollection collection = (GeometryCollection) geometry;
			for (int i = 0; i < collection.getNumGeometries(); i++) {
				addLinearGeometry(diagram, collection.getGeometryN(i), vertices);
			}
			return;
		}
		throw new IllegalArgumentException("Unsupported geometry type: " + geometry.getGeometryType());
	}

	private static void addPointSite(VoronoiDiagram diagram, Point point, List<Vertex> vertices) {
		Coordinate coordinate = point.getCoordinate();
		if (coordinate == null) {
			throw new IllegalArgumentException("Point geometry must contain a coordinate.");
		}
		vertices.add(diagram.insertPointSite(toPoint(coordinate)));
	}

	private static LineString normalizeClosedLineal(LineString lineString) {
		if (!lineString.isClosed()) {
			return lineString;
		}
		LineString normalized = (LineString) lineString.copy();
		normalized.normalize();
		return normalized;
	}

	private static Polygon normalizePolygon(Polygon polygon) {
		return (Polygon) polygon.norm();
	}
}
