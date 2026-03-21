package org.rogach.jopenvoronoi.jts;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.IdentityHashMap;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class JtsGeometryIOTest {

	private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();

	@Test
	public void polygonSitesSkipDuplicateClosingCoordinate() {
		Polygon polygon = GEOMETRY_FACTORY.createPolygon(new Coordinate[] { new Coordinate(-0.4, -0.4),
				new Coordinate(-0.4, 0.4), new Coordinate(0.4, 0.4), new Coordinate(0.4, -0.4),
				new Coordinate(-0.4, -0.4) });

		VoronoiDiagram diagram = new VoronoiDiagram();
		java.util.List<Vertex> vertices = JtsGeometryIO.addPolygonSites(diagram, polygon);

		assertEquals(4, vertices.size());
	}

	@Test
	public void multilineExportDeduplicatesTwinEdges() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		diagram.insertPointSite(-0.25, -0.1);
		diagram.insertPointSite(0.35, -0.05);
		diagram.insertPointSite(0.0, 0.4);

		MultiLineString exported = JtsGeometryIO.toMultiLineString(diagram);
		assertFalse(exported.isEmpty());

		IdentityHashMap<Edge, Boolean> uniqueBases = new IdentityHashMap<>();
		for (Edge edge : diagram.getDiagram().edges) {
			if (edge.valid && edge.base != null && edge.source != null && edge.target != null && edge.source.position != null
					&& edge.target.position != null) {
				uniqueBases.put(edge.base, Boolean.TRUE);
			}
		}
		assertEquals(uniqueBases.size(), exported.getNumGeometries());
	}

	@Test
	public void addGeometryNormalizesClosedLinealInputBeforeInsertion() {
		LineString ring = GEOMETRY_FACTORY.createLineString(new Coordinate[] { new Coordinate(1, 1), new Coordinate(1, 0),
				new Coordinate(0, 0), new Coordinate(0, 1), new Coordinate(1, 1) });
		LineString normalized = (LineString) ring.copy();
		normalized.normalize();

		VoronoiDiagram diagram = JtsGeometryIO.toVoronoiDiagram(ring);
		List<Vertex> vertices = collectInsertedVertices(diagram, normalized.getNumPoints() - 1);

		assertEquals(normalized.getNumPoints() - 1, vertices.size());
		for (int i = 0; i < vertices.size(); i++) {
			Coordinate expected = normalized.getCoordinateN(i);
			assertEquals(expected.getX(), vertices.get(i).position.x, 1e-12);
			assertEquals(expected.getY(), vertices.get(i).position.y, 1e-12);
		}
	}

	@Test
	public void addGeometryLeavesOpenLinealInputOrderUntouched() {
		LineString line = GEOMETRY_FACTORY
				.createLineString(new Coordinate[] { new Coordinate(0.1, 0.2), new Coordinate(0.3, 0.4), new Coordinate(0.5, 0.1) });

		VoronoiDiagram diagram = JtsGeometryIO.toVoronoiDiagram(line);
		List<Vertex> vertices = collectInsertedVertices(diagram, 3);

		assertEquals(3, vertices.size());
		for (int i = 0; i < vertices.size(); i++) {
			Coordinate expected = line.getCoordinateN(i);
			assertEquals(expected.getX(), vertices.get(i).position.x, 1e-12);
			assertEquals(expected.getY(), vertices.get(i).position.y, 1e-12);
		}
	}

	@Test
	public void addGeometryProcessesPointalMembersBeforeLinealMembers() {
		Point point = GEOMETRY_FACTORY.createPoint(new Coordinate(-0.7, -0.6));
		LineString line = GEOMETRY_FACTORY.createLineString(new Coordinate[] { new Coordinate(-0.2, -0.1), new Coordinate(0.2, 0.1) });
		Geometry mixed = GEOMETRY_FACTORY.createGeometryCollection(new Geometry[] { line, point });

		VoronoiDiagram diagram = JtsGeometryIO.toVoronoiDiagram(mixed);
		List<Vertex> vertices = collectInsertedVertices(diagram, 3);

		assertEquals(3, vertices.size());
		assertEquals(point.getX(), vertices.get(0).position.x, 1e-12);
		assertEquals(point.getY(), vertices.get(0).position.y, 1e-12);
	}

	private static List<Vertex> collectInsertedVertices(VoronoiDiagram diagram, int expectedCount) {
		List<Vertex> inserted = new java.util.ArrayList<>();
		for (Vertex vertex : diagram.getDiagram().vertices) {
			if (vertex.type == org.rogach.jopenvoronoi.vertex.VertexType.POINTSITE) {
				inserted.add(vertex);
			}
		}
		int from = Math.max(0, inserted.size() - expectedCount);
		return inserted.subList(from, inserted.size());
	}
}
