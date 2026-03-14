package org.rogach.jopenvoronoi;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.filter.MedialAxisFilter;
import org.rogach.jopenvoronoi.filter.PolygonInteriorFilter;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
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

	@Test
	public void medialAxisCanBeCollectedFromPolygonPointList() {
		VoronoiDiagram vd = new VoronoiDiagram();
		List<Point> polygon = List.of(new Point(-1.0, -1.0), new Point(1.0, -1.0), new Point(1.0, 1.0),
				new Point(-1.0, 1.0));

		List<Vertex> vertexHandles = new ArrayList<>();
		for (Point point : polygon) {
			vertexHandles.add(vd.insertPointSite(point));
		}
		for (int i = 0; i < vertexHandles.size(); i++) {
			vd.insertLineSite(vertexHandles.get(i), vertexHandles.get((i + 1) % vertexHandles.size()));
		}

		vd.filter(new PolygonInteriorFilter(true));
		vd.filter(new MedialAxisFilter());

		List<Edge> medialAxis = vd.getDiagram().edges.stream().filter(edge -> edge.valid)
				.filter(edge -> edge.type != EdgeType.LINESITE).filter(edge -> edge.type != EdgeType.NULLEDGE)
				.filter(edge -> edge.type != EdgeType.OUTEDGE).collect(Collectors.toList());

		Assertions.assertEquals(10, medialAxis.size(), "Expected the square to yield 10 medial-axis half-edges");
		Assertions.assertTrue(medialAxis.stream().allMatch(edge -> Math.abs(edge.source.position.x) <= 1.0
				&& Math.abs(edge.source.position.y) <= 1.0 && Math.abs(edge.target.position.x) <= 1.0
				&& Math.abs(edge.target.position.y) <= 1.0),
				"Medial-axis edges should stay inside the polygon bounds");
		Assertions.assertTrue(medialAxis.stream().anyMatch(edge -> edge.type == EdgeType.PARA_LINELINE),
				"Expected the filtered result to contain the square's center branch");
	}
}
