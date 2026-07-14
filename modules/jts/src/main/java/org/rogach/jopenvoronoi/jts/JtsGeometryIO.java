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
 * <p>
 * These are low-level conversion primitives. For the common derived products of
 * a polygon's Voronoi diagram — cells, interior cells, the dissolved medial
 * axis, and the medial-axis coverage — use {@link JtsVoronoi}.
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

		// The diagram requires every point site to be inserted before any line
		// site, so insertion is split into two global phases across the whole
		// geometry rather than done per element: first all point sites (standalone
		// points plus every linestring/ring vertex), then all line-site segments.
		// Doing this per element would fail for any input with more than one lineal
		// component (a MultiLineString, or a Polygon with holes).
		List<Vertex> pointSites = new ArrayList<>();
		addPointalGeometry(diagram, geometry, pointSites);

		List<LineString> lineals = new ArrayList<>();
		collectLinealComponents(geometry, lineals);

		List<List<Vertex>> linealVertices = new ArrayList<>(lineals.size());
		for (LineString lineal : lineals) {
			linealVertices.add(insertLineStringPoints(diagram, lineal));
		}
		for (int i = 0; i < lineals.size(); i++) {
			insertLineStringSegments(diagram, linealVertices.get(i), lineals.get(i).isClosed());
		}
	}

	public static List<Vertex> addLineStringSites(VoronoiDiagram diagram, LineString lineString) {
		if (diagram == null) {
			throw new IllegalArgumentException("Voronoi diagram cannot be null.");
		}
		List<Vertex> vertices = insertLineStringPoints(diagram, lineString);
		insertLineStringSegments(diagram, vertices, lineString.isClosed());
		return vertices;
	}

	public static List<Vertex> addPolygonSites(VoronoiDiagram diagram, Polygon polygon) {
		if (diagram == null) {
			throw new IllegalArgumentException("Voronoi diagram cannot be null.");
		}
		if (polygon == null) {
			throw new IllegalArgumentException("Polygon cannot be null.");
		}

		// Point sites for every ring must precede any line site, so insert all ring
		// vertices first and only then the segments (see addGeometry).
		List<LineString> rings = new ArrayList<>();
		rings.add(polygon.getExteriorRing());
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			rings.add(polygon.getInteriorRingN(i));
		}

		List<Vertex> vertices = new ArrayList<>();
		List<List<Vertex>> ringVertices = new ArrayList<>(rings.size());
		for (LineString ring : rings) {
			List<Vertex> inserted = insertLineStringPoints(diagram, ring);
			ringVertices.add(inserted);
			vertices.addAll(inserted);
		}
		for (int i = 0; i < rings.size(); i++) {
			insertLineStringSegments(diagram, ringVertices.get(i), rings.get(i).isClosed());
		}
		return vertices;
	}

	/**
	 * Inserts a point site for each distinct vertex of the linestring (dropping the
	 * duplicate closing coordinate of a closed ring) and returns them in order. No
	 * line sites are inserted, so this is safe to call for several linestrings
	 * before any segment is added.
	 */
	private static List<Vertex> insertLineStringPoints(VoronoiDiagram diagram, LineString lineString) {
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
			vertices.add(diagram.insertPointSite(toPoint(coordinates[i])));
		}
		return vertices;
	}

	/**
	 * Connects consecutive vertices with line sites, closing the loop when the
	 * originating linestring was closed. Assumes every vertex has already been
	 * inserted as a point site.
	 */
	private static void insertLineStringSegments(VoronoiDiagram diagram, List<Vertex> vertices, boolean closed) {
		for (int i = 0; i < vertices.size() - 1; i++) {
			diagram.insertLineSite(vertices.get(i), vertices.get(i + 1));
		}
		if (closed) {
			diagram.insertLineSite(vertices.get(vertices.size() - 1), vertices.get(0));
		}
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

	/**
	 * Collects the lineal components of a geometry — open/closed linestrings and
	 * polygon rings — as normalized linestrings, in insertion order, so that
	 * {@link #addGeometry} can insert all of their point sites before any line
	 * site. Standalone points are ignored (they are handled by
	 * {@link #addPointalGeometry}).
	 */
	private static void collectLinealComponents(Geometry geometry, List<LineString> out) {
		if (geometry instanceof Point) {
			return;
		}
		if (geometry instanceof Polygon) {
			Polygon polygon = normalizePolygon((Polygon) geometry);
			out.add(polygon.getExteriorRing());
			for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
				out.add(polygon.getInteriorRingN(i));
			}
			return;
		}
		if (geometry instanceof LineString) {
			out.add(normalizeClosedLineal((LineString) geometry));
			return;
		}
		if (geometry instanceof GeometryCollection) {
			GeometryCollection collection = (GeometryCollection) geometry;
			for (int i = 0; i < collection.getNumGeometries(); i++) {
				collectLinealComponents(collection.getGeometryN(i), out);
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
