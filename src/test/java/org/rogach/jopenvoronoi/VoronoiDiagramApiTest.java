package org.rogach.jopenvoronoi;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
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
		Assertions.assertEquals(diagram.numFaces(), diagram.getDiagram().numFaces());
		Assertions.assertFalse(diagram.getDiagram().faceEdges(diagram.getFaces().get(0)).isEmpty());
	}
}
