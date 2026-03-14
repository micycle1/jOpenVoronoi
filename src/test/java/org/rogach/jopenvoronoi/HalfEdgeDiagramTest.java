package org.rogach.jopenvoronoi;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.util.Pair;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class HalfEdgeDiagramTest {

	@Test
	public void addTwinEdgesUsesFirstHalfEdgeAsCanonicalBase() {
		HalfEdgeDiagram diagram = new HalfEdgeDiagram();
		Vertex source = diagram.add_vertex();
		Vertex target = diagram.add_vertex();

		Pair<Edge, Edge> pair = diagram.add_twin_edges(source, target);
		Edge edge = pair.getFirst();
		Edge twin = pair.getSecond();

		Assertions.assertSame(twin, edge.twin);
		Assertions.assertSame(edge, twin.twin);
		Assertions.assertSame(edge, edge.base);
		Assertions.assertSame(edge, twin.base);
	}
}
