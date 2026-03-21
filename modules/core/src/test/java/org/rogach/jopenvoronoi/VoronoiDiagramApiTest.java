package org.rogach.jopenvoronoi;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
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
