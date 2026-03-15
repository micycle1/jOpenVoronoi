package org.rogach.jopenvoronoi;

import java.util.Map.Entry;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.vertex.Vertex;
import org.rogach.jopenvoronoi.vertex.VertexStatus;
import org.rogach.jopenvoronoi.vertex.VertexType;

public class EdgeTest {

	@Test
	public void pointAndRadiusInterpolatesLinearEdges() {
		Edge edge = new Edge(new Vertex(new Point(0, 0), VertexStatus.UNDECIDED, VertexType.NORMAL, 0.0),
				new Vertex(new Point(2, 0), VertexStatus.UNDECIDED, VertexType.NORMAL, 0.0));
		edge.type = EdgeType.LINE;

		Face face = new Face();
		face.setSite(new PointSite(new Point(1, 1), face));
		edge.face = face;

		Entry<Point, Double> sample = edge.pointAndRadius(0.25);

		assertPointEquals(new Point(0.5, 0.0), sample.getKey());
		Assertions.assertEquals(Math.sqrt(1.25), sample.getValue(), 1e-12);
	}

	@Test
	@SuppressWarnings("deprecation")
	public void edgePointDelegatesToPointAndRadiusForParabolas() {
		PointSite pointSite = new PointSite(new Point(0, 0));
		LineSite lineSite = new LineSite(new Point(-2, 1), new Point(2, 1), 1.0);

		Edge edge = new Edge(new Vertex(), new Vertex());
		edge.setParameters(pointSite, lineSite, true);

		double startRadius = 1.0;
		double endRadius = 3.0;
		edge.source = new Vertex(edge.point(startRadius), VertexStatus.UNDECIDED, VertexType.NORMAL, startRadius);
		edge.target = new Vertex(edge.point(endRadius), VertexStatus.UNDECIDED, VertexType.NORMAL, endRadius);

		Entry<Point, Double> sample = edge.pointAndRadius(0.5);
		Entry<Point, Double> deprecatedSample = edge.edgePoint(0.5);

		assertPointEquals(edge.point(2.0), sample.getKey());
		Assertions.assertEquals(2.0, sample.getValue(), 1e-12);
		assertPointEquals(sample.getKey(), deprecatedSample.getKey());
		Assertions.assertEquals(sample.getValue(), deprecatedSample.getValue(), 1e-12);
	}

	private static void assertPointEquals(Point expected, Point actual) {
		Assertions.assertEquals(expected.x, actual.x, 1e-12);
		Assertions.assertEquals(expected.y, actual.y, 1e-12);
	}
}
