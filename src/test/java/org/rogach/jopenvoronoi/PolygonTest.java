package org.rogach.jopenvoronoi;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class PolygonTest {

	@Test
	public void polygonTest() {
		// create diagram
		VoronoiDiagram vd = new VoronoiDiagram();

		// simple square polygon
		double a = 0.4567;

		// keep handles returned by inserting point-sites so we can create line-sites
		List<Vertex> vertexHandles = new ArrayList<>();
		vertexHandles.add(vd.insertPointSite(new Point(-a, -a)));
		vertexHandles.add(vd.insertPointSite(new Point(-a, a)));
		vertexHandles.add(vd.insertPointSite(new Point(a, a)));
		vertexHandles.add(vd.insertPointSite(new Point(a, -a)));

		// insert line-sites connecting the points into a closed polygon
		for (int i = 0; i < vertexHandles.size(); i++) {
			int next = (i == vertexHandles.size() - 1) ? 0 : i + 1;
			vd.insertLineSite(vertexHandles.get(i), vertexHandles.get(next));
		}

		// basic checks
		Assertions.assertTrue(vd.check(), "VoronoiDiagram.check() should return true");

		// sanity on site counts
		Assertions.assertEquals(4, vd.numPointSites(), "Should have 4 point sites");
		Assertions.assertEquals(4, vd.numLineSites(), "Should have 4 line sites");

		// retrieve diagram and ensure it's present
		HalfEdgeDiagram diagram = vd.getDiagram();
		Assertions.assertNotNull(diagram, "HalfEdgeDiagram should not be null");

		// Optionally iterate faces/edges to exercise the API
		vd.getFaces().forEach(face -> {
			Point pos = face.site.position();
			// ensure site position is present (just a basic non-null check)
			Assertions.assertNotNull(pos);
			diagram.faceEdges(face).forEach(edge -> {
				Assertions.assertNotNull(edge.source);
				Assertions.assertNotNull(edge.target);
				// ensure source/target positions exist
				Assertions.assertNotNull(edge.source.position);
				Assertions.assertNotNull(edge.target.position);
			});
		});
	}
}
