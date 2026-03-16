package org.rogach.jopenvoronoi.pocket;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.filter.MedialAxisFilter;
import org.rogach.jopenvoronoi.filter.PolygonInteriorFilter;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class PocketPathTest {

	@Test
	public void emptyMicListProducesEmptyPath() {
		List<PocketPath.Segment> path = PocketPath.toPath(new ArrayList<>());
		Assertions.assertTrue(path.isEmpty());
	}

	@Test
	public void nullMicListProducesEmptyPath() {
		List<PocketPath.Segment> path = PocketPath.toPath(null);
		Assertions.assertTrue(path.isEmpty());
	}

	@Test
	public void singleMicProducesEmptyPath() {
		// A single MIC is just the starting circle — no toolpath segments yet
		MIC mic = new MIC();
		mic.c1 = new Point(0, 0);
		mic.c2 = new Point(0, 0);
		mic.r1 = 1.0;
		mic.r2 = 1.0;
		List<MIC> mics = new ArrayList<>();
		mics.add(mic);
		List<PocketPath.Segment> path = PocketPath.toPath(mics);
		Assertions.assertTrue(path.isEmpty());
	}

	@Test
	public void twoMicsProducesFourSegments() {
		// Build two MICs with known tangent points
		MIC mic0 = new MIC();
		mic0.c2 = new Point(0, 0);
		mic0.r2 = 1.0;

		MIC mic1 = new MIC();
		mic1.c1 = new Point(0, 0);
		mic1.c2 = new Point(2, 0);
		mic1.r1 = 1.0;
		mic1.r2 = 1.0;
		// With equal radii along x-axis, tangent points at y = ±1
		mic1.t1 = new Point(0, -1);
		mic1.t2 = new Point(0, 1);
		mic1.t3 = new Point(2, -1);
		mic1.t4 = new Point(2, 1);
		mic1.newBranch = false;

		List<MIC> mics = new ArrayList<>();
		mics.add(mic0);
		mics.add(mic1);

		List<PocketPath.Segment> path = PocketPath.toPath(mics);
		// Should produce: arc on c1, line t1→t3, arc on c2, line t4→t2
		Assertions.assertEquals(4, path.size());

		// First segment: arc on previous circle
		Assertions.assertTrue(path.get(0).isArc());
		Assertions.assertEquals(1.0, path.get(0).radius, 1e-12);

		// Second segment: line from t1 to t3
		Assertions.assertTrue(path.get(1).isLine());
		assertPointEquals(new Point(0, -1), path.get(1).start);
		assertPointEquals(new Point(2, -1), path.get(1).end);

		// Third segment: arc on new circle
		Assertions.assertTrue(path.get(2).isArc());
		Assertions.assertEquals(1.0, path.get(2).radius, 1e-12);

		// Fourth segment: line from t4 to t2
		Assertions.assertTrue(path.get(3).isLine());
		assertPointEquals(new Point(2, 1), path.get(3).start);
		assertPointEquals(new Point(0, 1), path.get(3).end);
	}

	@Test
	public void branchMicInsertsRapid() {
		MIC mic0 = new MIC();
		mic0.c2 = new Point(0, 0);
		mic0.r2 = 1.0;

		MIC mic1 = new MIC();
		mic1.c1 = new Point(0, 0);
		mic1.c2 = new Point(2, 0);
		mic1.r1 = 1.0;
		mic1.r2 = 0.5;
		mic1.t1 = new Point(0, -1);
		mic1.t2 = new Point(0, 1);
		mic1.t3 = new Point(2, -0.5);
		mic1.t4 = new Point(2, 0.5);
		mic1.newBranch = true;
		mic1.cPrev = new Point(-1, -1);
		mic1.rPrev = 0.8;

		List<MIC> mics = new ArrayList<>();
		mics.add(mic0);
		mics.add(mic1);

		List<PocketPath.Segment> path = PocketPath.toPath(mics);
		// First segment should be a rapid move
		Assertions.assertTrue(path.get(0).isRapid());
		assertPointEquals(new Point(-1, -1), path.get(0).end);
	}

	@Test
	public void squarePocketProducesNonEmptyPath() {
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

		for (List<MIC> comp : pocket.getMicComponents()) {
			List<PocketPath.Segment> path = PocketPath.toPath(comp);
			Assertions.assertFalse(path.isEmpty(), "toolpath should not be empty");

			// Verify all segments have valid coordinates
			for (PocketPath.Segment seg : path) {
				Assertions.assertNotNull(seg.start, "segment start must not be null");
				Assertions.assertNotNull(seg.end, "segment end must not be null");
				Assertions.assertFalse(Double.isNaN(seg.start.x));
				Assertions.assertFalse(Double.isNaN(seg.start.y));
				Assertions.assertFalse(Double.isNaN(seg.end.x));
				Assertions.assertFalse(Double.isNaN(seg.end.y));
				if (seg.isArc()) {
					Assertions.assertNotNull(seg.center, "arc center must not be null");
					Assertions.assertTrue(seg.radius > 0, "arc radius must be positive");
				}
			}
		}
	}

	@Test
	public void segmentTypesAreCorrect() {
		// Build a simple two-MIC chain to get all segment types
		MIC mic0 = new MIC();
		mic0.c2 = new Point(0, 0);
		mic0.r2 = 1.0;

		MIC mic1 = new MIC();
		mic1.c1 = new Point(0, 0);
		mic1.c2 = new Point(2, 0);
		mic1.r1 = 1.0;
		mic1.r2 = 1.0;
		mic1.t1 = new Point(0, -1);
		mic1.t2 = new Point(0, 1);
		mic1.t3 = new Point(2, -1);
		mic1.t4 = new Point(2, 1);
		mic1.newBranch = true;
		mic1.cPrev = new Point(-1, -1);
		mic1.rPrev = 0.8;

		List<MIC> mics = new ArrayList<>();
		mics.add(mic0);
		mics.add(mic1);

		List<PocketPath.Segment> path = PocketPath.toPath(mics);
		// First segment is a rapid (newBranch), followed by arc, line, arc, line
		Assertions.assertTrue(path.get(0).isRapid());
		Assertions.assertFalse(path.get(0).isLine());
		Assertions.assertFalse(path.get(0).isArc());
		Assertions.assertTrue(path.get(1).isArc());
		Assertions.assertTrue(path.get(2).isLine());
		Assertions.assertTrue(path.get(3).isArc());
		Assertions.assertTrue(path.get(4).isLine());

		// Verify toString doesn't throw
		for (PocketPath.Segment seg : path) {
			Assertions.assertNotNull(seg.toString());
		}
	}

	private static void assertPointEquals(Point expected, Point actual) {
		Assertions.assertEquals(expected.x, actual.x, 1e-12);
		Assertions.assertEquals(expected.y, actual.y, 1e-12);
	}
}
