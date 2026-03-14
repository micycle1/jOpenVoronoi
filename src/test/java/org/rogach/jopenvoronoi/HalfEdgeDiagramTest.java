package org.rogach.jopenvoronoi;

import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.util.Pair;
import org.rogach.jopenvoronoi.vertex.Vertex;
import org.rogach.jopenvoronoi.vertex.VertexStatus;
import org.rogach.jopenvoronoi.vertex.VertexType;

public class HalfEdgeDiagramTest {

	@Test
	public void setNextCycleMaintainsPrevLinks() {
		HalfEdgeDiagram diagram = new HalfEdgeDiagram();
		Face face = diagram.add_face();
		Vertex a = diagram.add_vertex(vertex(0.0, 0.0, VertexType.NORMAL));
		Vertex b = diagram.add_vertex(vertex(1.0, 0.0, VertexType.NORMAL));
		Vertex c = diagram.add_vertex(vertex(0.0, 1.0, VertexType.NORMAL));
		Edge ab = diagram.add_edge(a, b);
		Edge bc = diagram.add_edge(b, c);
		Edge ca = diagram.add_edge(c, a);

		diagram.set_next_cycle(List.of(ab, bc, ca), face, 1.0);

		assertPrevLinks(diagram, face);
		Assertions.assertSame(ca, diagram.previous_edge(ab));
		Assertions.assertSame(ab, diagram.previous_edge(bc));
		Assertions.assertSame(bc, diagram.previous_edge(ca));
	}

	@Test
	public void splitAndRemoveEdgeKeepsPrevLinksConsistent() {
		DiagramFixture fixture = createTwoFaceFixture();
		HalfEdgeDiagram diagram = fixture.diagram;

		Vertex split = diagram.add_vertex(vertex(1.0, 0.0, VertexType.SPLIT));
		diagram.add_vertex_in_edge(split, fixture.shared);

		Assertions.assertEquals(8, diagram.num_edges());
		Assertions.assertFalse(diagram.edges.contains(fixture.shared));
		Assertions.assertFalse(diagram.edges.contains(fixture.sharedTwin));
		assertPrevLinks(diagram, fixture.leftFace);
		assertPrevLinks(diagram, fixture.rightFace);

		Edge aToSplit = diagram.edge(fixture.a, split);
		Edge splitToB = diagram.edge(split, fixture.b);
		Edge bToSplit = diagram.edge(fixture.b, split);
		Edge splitToA = diagram.edge(split, fixture.a);
		Assertions.assertSame(fixture.leftFace.edge, aToSplit);
		Assertions.assertSame(fixture.rightFace.edge, bToSplit);
		Assertions.assertSame(aToSplit, splitToB.prev);
		Assertions.assertSame(bToSplit, splitToA.prev);

		diagram.remove_deg2_vertex(split);

		Assertions.assertEquals(6, diagram.num_edges());
		Assertions.assertFalse(diagram.vertices.contains(split));
		Assertions.assertTrue(diagram.has_edge(fixture.a, fixture.b));
		Assertions.assertTrue(diagram.has_edge(fixture.b, fixture.a));
		assertPrevLinks(diagram, fixture.leftFace);
		assertPrevLinks(diagram, fixture.rightFace);
	}

	private static void assertPrevLinks(HalfEdgeDiagram diagram, Face face) {
		List<Edge> edges = diagram.faceEdges(face);
		for (int i = 0; i < edges.size(); i++) {
			Edge edge = edges.get(i);
			Edge expectedPrev = edges.get((i + edges.size() - 1) % edges.size());
			Assertions.assertSame(expectedPrev, edge.prev);
			Assertions.assertSame(expectedPrev, diagram.previous_edge(edge));
		}
	}

	private static DiagramFixture createTwoFaceFixture() {
		HalfEdgeDiagram diagram = new HalfEdgeDiagram();
		Face leftFace = diagram.add_face();
		Face rightFace = diagram.add_face();
		Vertex a = diagram.add_vertex(vertex(0.0, 0.0, VertexType.NORMAL));
		Vertex b = diagram.add_vertex(vertex(2.0, 0.0, VertexType.NORMAL));
		Vertex top = diagram.add_vertex(vertex(1.0, 1.0, VertexType.NORMAL));
		Vertex bottom = diagram.add_vertex(vertex(1.0, -1.0, VertexType.NORMAL));

		Edge topToA = diagram.add_edge(top, a);
		Pair<Edge, Edge> sharedEdges = diagram.add_twin_edges(a, b);
		Edge bToTop = diagram.add_edge(b, top);
		diagram.set_next_cycle(List.of(topToA, sharedEdges.getFirst(), bToTop), leftFace, 1.0);

		Edge bottomToB = diagram.add_edge(bottom, b);
		Edge aToBottom = diagram.add_edge(a, bottom);
		diagram.set_next_cycle(List.of(bottomToB, sharedEdges.getSecond(), aToBottom), rightFace, -1.0);

		return new DiagramFixture(diagram, leftFace, rightFace, a, b, sharedEdges.getFirst(), sharedEdges.getSecond());
	}

	private static Vertex vertex(double x, double y, VertexType type) {
		return new Vertex(new Point(x, y), VertexStatus.OUT, type);
	}

	private static class DiagramFixture {
		private final HalfEdgeDiagram diagram;
		private final Face leftFace;
		private final Face rightFace;
		private final Vertex a;
		private final Vertex b;
		private final Edge shared;
		private final Edge sharedTwin;

		private DiagramFixture(HalfEdgeDiagram diagram, Face leftFace, Face rightFace, Vertex a, Vertex b, Edge shared,
				Edge sharedTwin) {
			this.diagram = diagram;
			this.leftFace = leftFace;
			this.rightFace = rightFace;
			this.a = a;
			this.b = b;
			this.shared = shared;
			this.sharedTwin = sharedTwin;
		}
	}
}
