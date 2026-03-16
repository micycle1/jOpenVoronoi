package org.rogach.jopenvoronoi;

import java.util.List;
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
	public void micSampleInterpolatesLinearEdges() {
		Edge edge = new Edge(new Vertex(new Point(0, 0), VertexStatus.UNDECIDED, VertexType.NORMAL, 0.0),
				new Vertex(new Point(2, 0), VertexStatus.UNDECIDED, VertexType.NORMAL, 0.0));
		edge.type = EdgeType.LINE;

		Face face = new Face();
		face.setSite(new PointSite(new Point(1, 1), face));
		edge.face = face;

		Entry<Point, Double> sample = edge.micSample(0.25);

		assertPointEquals(new Point(0.5, 0.0), sample.getKey());
		Assertions.assertEquals(Math.sqrt(1.25), sample.getValue(), 1e-12);
	}

	@Test
	public void micSampleInterpolatesParabolaRadius() {
		PointSite pointSite = new PointSite(new Point(0, 0));
		LineSite lineSite = new LineSite(new Point(-2, 1), new Point(2, 1), 1.0);

		Edge edge = new Edge(new Vertex(), new Vertex());
		edge.setParameters(pointSite, lineSite, true);

		double startRadius = 1.0;
		double endRadius = 3.0;
		edge.source = new Vertex(edge.point(startRadius), VertexStatus.UNDECIDED, VertexType.NORMAL, startRadius);
		edge.target = new Vertex(edge.point(endRadius), VertexStatus.UNDECIDED, VertexType.NORMAL, endRadius);

		Entry<Point, Double> sample = edge.micSample(0.5);

		assertPointEquals(edge.point(2.0), sample.getKey());
		Assertions.assertEquals(2.0, sample.getValue(), 1e-12);
	}

	private static void assertPointEquals(Point expected, Point actual) {
		Assertions.assertEquals(expected.x, actual.x, 1e-12);
		Assertions.assertEquals(expected.y, actual.y, 1e-12);
	}

	@Test
	public void samplePointsLinearEdgeIncludesStartAndEnd() {
		// PP bisector of (0,0) and (2,0) is the vertical line x=1.
		// point(t) for t >= 1 gives (1, sqrt(t^2-1)) [sign=true].
		PointSite s1 = new PointSite(new Point(0, 0));
		PointSite s2 = new PointSite(new Point(2, 0));
		Edge edge = new Edge(new Vertex(), new Vertex());
		edge.setParameters(s1, s2, true);

		double startT = 1.0;              // apex
		double endT = Math.sqrt(2.0);     // t giving (1, 1)
		edge.source = new Vertex(edge.point(startT), VertexStatus.UNDECIDED, VertexType.NORMAL, startT);
		edge.target = new Vertex(edge.point(endT), VertexStatus.UNDECIDED, VertexType.NORMAL, endT);

		List<Point> pts = edge.samplePoints(3);

		Assertions.assertEquals(3, pts.size());
		assertPointEquals(edge.point(startT), pts.get(0));
		assertPointEquals(edge.point((startT + endT) / 2), pts.get(1));
		assertPointEquals(edge.point(endT), pts.get(2));
	}

	@Test
	public void samplePointsParabolaEdgeIncludesStartAndEnd() {
		PointSite pointSite = new PointSite(new Point(0, 0));
		LineSite lineSite = new LineSite(new Point(-2, 1), new Point(2, 1), 1.0);

		Edge edge = new Edge(new Vertex(), new Vertex());
		edge.setParameters(pointSite, lineSite, true);

		double startRadius = 1.0;
		double endRadius = 3.0;
		edge.source = new Vertex(edge.point(startRadius), VertexStatus.UNDECIDED, VertexType.NORMAL, startRadius);
		edge.target = new Vertex(edge.point(endRadius), VertexStatus.UNDECIDED, VertexType.NORMAL, endRadius);

		List<Point> pts = edge.samplePoints(3);

		Assertions.assertEquals(3, pts.size());
		assertPointEquals(edge.point(startRadius), pts.get(0));
		assertPointEquals(edge.point(2.0), pts.get(1));
		assertPointEquals(edge.point(endRadius), pts.get(2));
	}

	@Test
	public void samplePointsReturnsExactlyTwoPointsWhenNIsTwo() {
		PointSite s1 = new PointSite(new Point(0, 0));
		PointSite s2 = new PointSite(new Point(2, 0));
		Edge edge = new Edge(new Vertex(), new Vertex());
		edge.setParameters(s1, s2, true);

		double startT = 1.0;
		double endT = Math.sqrt(2.0);
		edge.source = new Vertex(edge.point(startT), VertexStatus.UNDECIDED, VertexType.NORMAL, startT);
		edge.target = new Vertex(edge.point(endT), VertexStatus.UNDECIDED, VertexType.NORMAL, endT);

		List<Point> pts = edge.samplePoints(2);

		Assertions.assertEquals(2, pts.size());
		assertPointEquals(edge.point(startT), pts.get(0));
		assertPointEquals(edge.point(endT), pts.get(1));
	}

	@Test
	public void samplePointsThrowsForNLessThanTwo() {
		PointSite s1 = new PointSite(new Point(0, 0));
		PointSite s2 = new PointSite(new Point(2, 0));
		Edge edge = new Edge(new Vertex(), new Vertex());
		edge.setParameters(s1, s2, true);
		edge.source = new Vertex(edge.point(1.0), VertexStatus.UNDECIDED, VertexType.NORMAL, 1.0);
		edge.target = new Vertex(edge.point(Math.sqrt(2.0)), VertexStatus.UNDECIDED, VertexType.NORMAL, Math.sqrt(2.0));

		Assertions.assertThrows(IllegalArgumentException.class, () -> edge.samplePoints(1));
		Assertions.assertThrows(IllegalArgumentException.class, () -> edge.samplePoints(0));
	}
}
