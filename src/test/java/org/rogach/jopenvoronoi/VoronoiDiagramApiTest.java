package org.rogach.jopenvoronoi;

import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.util.Pair;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class VoronoiDiagramApiTest {

	@Test
	public void insertPointSiteRejectsOutOfBoundsPoints() {
		VoronoiDiagram diagram = new VoronoiDiagram(1.0);

		IllegalArgumentException error = Assertions.assertThrows(IllegalArgumentException.class,
				() -> diagram.insertPointSite(new Point(1.0, 0.0)));

		Assertions.assertTrue(error.getMessage().contains("far radius"));
	}

	@Test
	public void insertPointSiteRejectsDuplicates() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		diagram.insertPointSite(new Point(0.25, 0.25));

		IllegalArgumentException error = Assertions.assertThrows(IllegalArgumentException.class,
				() -> diagram.insertPointSite(new Point(0.25, 0.25)));

		Assertions.assertTrue(error.getMessage().contains("duplicate"));
	}

	@Test
	public void insertPointSiteRejectsInsertionsAfterLineSites() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		Vertex start = diagram.insertPointSite(new Point(-0.2, 0.0));
		Vertex end = diagram.insertPointSite(new Point(0.2, 0.0));
		diagram.insertLineSite(start, end);

		IllegalStateException error = Assertions.assertThrows(IllegalStateException.class,
				() -> diagram.insertPointSite(new Point(0.0, 0.3)));

		Assertions.assertTrue(error.getMessage().contains("before inserting any line sites"));
	}

	@Test
	public void insertPointSitesReturnsHandlesInIterationOrder() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		List<Point> points = List.of(new Point(-0.4, -0.4), new Point(-0.4, 0.4), new Point(0.4, 0.4));

		List<Vertex> handles = diagram.insertPointSites(points);

		Assertions.assertEquals(points.size(), handles.size());
		for (int i = 0; i < points.size(); i++) {
			Assertions.assertEquals(points.get(i), handles.get(i).position);
		}
		Assertions.assertEquals(points.size(), diagram.numPointSites());
	}

	@Test
	public void insertLineSegmentsReusesExistingAndSharedEndpoints() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		Point first = new Point(-0.4, -0.4);
		Point second = new Point(-0.4, 0.4);
		Point third = new Point(0.4, 0.4);
		Point fourth = new Point(0.4, -0.4);
		diagram.insertPointSite(first);

		diagram.insertLineSegments(List.of(new Pair<>(first, second), new Pair<>(second, third), new Pair<>(third, fourth),
				new Pair<>(fourth, first)));

		Assertions.assertEquals(4, diagram.numPointSites());
		Assertions.assertEquals(4, diagram.numLineSites());
		Assertions.assertTrue(diagram.check(), "VoronoiDiagram.check() should return true after bulk segment insertion");
	}

	@Test
	public void camelCaseApiRemainsUsable() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		Vertex first = diagram.insertPointSite(new Point(-0.4, -0.4));
		Vertex second = diagram.insertPointSite(new Point(-0.4, 0.4));
		Vertex third = diagram.insertPointSite(new Point(0.4, 0.4));
		Vertex fourth = diagram.insertPointSite(new Point(0.4, -0.4));

		diagram.insertLineSite(first, second);
		diagram.insertLineSite(second, third);
		diagram.insertLineSite(third, fourth);
		diagram.insertLineSite(fourth, first);

		Assertions.assertEquals(4, diagram.numPointSites());
		Assertions.assertEquals(4, diagram.numLineSites());
		Assertions.assertEquals(diagram.getFaces().size(), diagram.numFaces());
		Assertions.assertTrue(diagram.numAllFaces() >= diagram.numFaces());
		Assertions.assertTrue(diagram.getAllFaces().size() >= diagram.getFaces().size());
		Assertions.assertSame(diagram.getFace(first), first.face);

		Face face = diagram.getFaces().stream().filter(Face::isPointSiteFace).findFirst().orElseThrow();
		Assertions.assertFalse(face.getEdges().isEmpty());
		Assertions.assertFalse(face.getVertices().isEmpty());
		Assertions.assertEquals(face.getEdges(), diagram.getFaceEdges(face));
		Assertions.assertEquals(face.getVertices(), diagram.getFaceVertices(face));
		Assertions.assertFalse(diagram.getAdjacentFaces(face).isEmpty());
		Assertions.assertTrue(diagram.getFaces().stream().noneMatch(Face::isNullFace));
		Assertions.assertEquals(diagram.getNonNullFaces().size(), diagram.getFaces().size());
		Assertions.assertTrue(diagram.getAllFaces().stream().anyMatch(Face::isNullFace));
		Assertions.assertTrue(diagram.check());
	}
}
