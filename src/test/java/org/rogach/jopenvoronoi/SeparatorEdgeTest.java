package org.rogach.jopenvoronoi;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class SeparatorEdgeTest {

	@Test
	public void separatorEdgesKeepPointSideKPositive() {
		VoronoiDiagram diagram = new VoronoiDiagram();
		Vertex start = diagram.insertPointSite(new Point(-0.2, 0.0));
		Vertex end = diagram.insertPointSite(new Point(0.2, 0.0));

		diagram.insertLineSite(start, end);

		List<Edge> separatorEdges = diagram.getDiagram().edges.stream().filter(edge -> edge.type == EdgeType.SEPARATOR)
				.collect(Collectors.toList());
		assertFalse(separatorEdges.isEmpty());

		for (Edge edge : separatorEdges) {
			if (edge.face != null && edge.face.site != null && edge.face.site.isPoint()) {
				assertEquals(1.0, edge.k, 0.0);
			}
			if (edge.twin.face != null && edge.twin.face.site != null && edge.twin.face.site.isPoint()) {
				assertEquals(1.0, edge.twin.k, 0.0);
			}
		}
	}
}
