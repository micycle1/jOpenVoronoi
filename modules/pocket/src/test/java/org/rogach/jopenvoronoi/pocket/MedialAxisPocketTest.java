package org.rogach.jopenvoronoi.pocket;

import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.filter.MedialAxisFilter;
import org.rogach.jopenvoronoi.filter.PolygonInteriorFilter;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class MedialAxisPocketTest {

	@Test
	public void squarePocketProducesMics() {
		VoronoiDiagram vd = new VoronoiDiagram();
		Vertex v1 = vd.insertPointSite(new Point(-0.4, -0.4));
		Vertex v2 = vd.insertPointSite(new Point(-0.4, 0.4));
		Vertex v3 = vd.insertPointSite(new Point(0.4, 0.4));
		Vertex v4 = vd.insertPointSite(new Point(0.4, -0.4));
		vd.insertLineSite(v1, v2);
		vd.insertLineSite(v2, v3);
		vd.insertLineSite(v3, v4);
		vd.insertLineSite(v4, v1);

		vd.filter(new PolygonInteriorFilter(true));
		vd.filter(new MedialAxisFilter());

		MedialAxisPocket pocket = new MedialAxisPocket(vd.getDiagram());
		pocket.setWidth(0.1);
		pocket.run();

		List<List<MIC>> components = pocket.getMicComponents();
		Assertions.assertFalse(components.isEmpty(), "should produce at least one component");

		int totalMics = components.stream().mapToInt(List::size).sum();
		Assertions.assertTrue(totalMics > 0, "should produce at least one MIC");

		// verify all MIC radii are positive and centers have valid coordinates
		for (List<MIC> comp : components) {
			for (MIC mic : comp) {
				Assertions.assertTrue(mic.r2 > 0, "MIC radius should be positive");
				Assertions.assertNotNull(mic.c2, "MIC center should not be null");
				Assertions.assertFalse(Double.isNaN(mic.c2.x), "MIC center x should not be NaN");
				Assertions.assertFalse(Double.isNaN(mic.c2.y), "MIC center y should not be NaN");
			}
		}
	}

	@Test
	public void cutWidthComputation() {
		Point c1 = new Point(0, 0);
		Point c2 = new Point(1, 0);
		double w = MedialAxisPocket.cutWidth(c1, 0.5, c2, 0.7);
		// |c2 - c1| + r2 - r1 = 1.0 + 0.7 - 0.5 = 1.2
		Assertions.assertEquals(1.2, w, 1e-12);
	}

	@Test
	public void bitangentPointsEqualRadii() {
		Point c1 = new Point(0, 0);
		Point c2 = new Point(2, 0);
		List<Point> pts = MedialAxisPocket.bitangentPoints(c1, 1.0, c2, 1.0);
		Assertions.assertEquals(4, pts.size());

		// With equal radii along x-axis, tangent points should be at y = ±1
		Assertions.assertEquals(0.0, pts.get(0).x, 1e-12);
		Assertions.assertEquals(-1.0, pts.get(0).y, 1e-12);
		Assertions.assertEquals(0.0, pts.get(1).x, 1e-12);
		Assertions.assertEquals(1.0, pts.get(1).y, 1e-12);
		Assertions.assertEquals(2.0, pts.get(2).x, 1e-12);
		Assertions.assertEquals(-1.0, pts.get(2).y, 1e-12);
		Assertions.assertEquals(2.0, pts.get(3).x, 1e-12);
		Assertions.assertEquals(1.0, pts.get(3).y, 1e-12);
	}

	@Test
	public void trianglePocketProducesMics() {
		VoronoiDiagram vd = new VoronoiDiagram();
		Vertex v1 = vd.insertPointSite(new Point(0.0, 0.4));
		Vertex v2 = vd.insertPointSite(new Point(-0.4, -0.3));
		Vertex v3 = vd.insertPointSite(new Point(0.4, -0.3));
		vd.insertLineSite(v1, v2);
		vd.insertLineSite(v2, v3);
		vd.insertLineSite(v3, v1);

		vd.filter(new PolygonInteriorFilter(true));
		vd.filter(new MedialAxisFilter());

		MedialAxisPocket pocket = new MedialAxisPocket(vd.getDiagram());
		pocket.setWidth(0.05);
		pocket.run();

		List<List<MIC>> components = pocket.getMicComponents();
		Assertions.assertFalse(components.isEmpty());

		int totalMics = components.stream().mapToInt(List::size).sum();
		Assertions.assertTrue(totalMics > 0);
	}

	@Test
	public void emptyDiagramProducesNoComponents() {
		// An empty HalfEdgeDiagram with no edges should produce no components
		org.rogach.jopenvoronoi.HalfEdgeDiagram emptyGraph = new org.rogach.jopenvoronoi.HalfEdgeDiagram();
		MedialAxisPocket pocket = new MedialAxisPocket(emptyGraph);
		pocket.run();

		Assertions.assertTrue(pocket.getMicComponents().isEmpty());
	}
}
