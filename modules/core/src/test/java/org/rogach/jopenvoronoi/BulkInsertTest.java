package org.rogach.jopenvoronoi;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class BulkInsertTest {

	private static List<Point> randomPoints(int n, long seed) {
		var rnd = new Random(seed);
		var pts = new ArrayList<Point>(n);
		for (var i = 0; i < n; i++) {
			pts.add(new Point(rnd.nextDouble() * 1000 - 500, rnd.nextDouble() * 1000 - 500));
		}
		return pts;
	}

	@Test
	public void bulkMatchesSingleInsertion() {
		var pts = randomPoints(5000, 42);

		var single = new VoronoiDiagram();
		for (Point p : pts) {
			single.insertPointSite(p);
		}
		var bulk = new VoronoiDiagram();
		bulk.insertPointSites(pts);

		assertEquals(single.numPointSites(), bulk.numPointSites());
		assertEquals(single.numVertices(), bulk.numVertices());
		assertEquals(single.numFaces(), bulk.numFaces());
		assertEquals(single.numAllFaces(), bulk.numAllFaces());
		assertTrue(single.check());
		assertTrue(bulk.check());
	}

	@Test
	public void resultIsAlignedToInputOrder() {
		var pts = randomPoints(500, 7);
		var vd = new VoronoiDiagram();
		var handles = vd.insertPointSites(pts);
		assertEquals(pts.size(), handles.size());
		for (var i = 0; i < pts.size(); i++) {
			assertEquals(pts.get(i).x, handles.get(i).position.x);
			assertEquals(pts.get(i).y, handles.get(i).position.y);
		}
	}

	@Test
	public void duplicatesResolveToSameVertex() {
		var vd = new VoronoiDiagram();
		var p = new Point(10, 20);
		var handles = vd.insertPointSites(List.of(p, new Point(30, 40), new Point(10, 20)));
		assertSame(handles.get(0), handles.get(2));
		assertEquals(2, vd.numPointSites());
	}

	@Test
	public void emptyInputIsANoOp() {
		var vd = new VoronoiDiagram();
		assertTrue(vd.insertPointSites(new ArrayList<>()).isEmpty());
		assertEquals(0, vd.numPointSites());
	}

	@Test
	public void failFastValidationLeavesDiagramUnmutated() {
		var vd = new VoronoiDiagram();
		// invalid point at the END of the list: nothing must be inserted
		var bad = new ArrayList<Point>();
		bad.add(new Point(1, 1));
		bad.add(new Point(99999, 0)); // outside farRadius 5000
		assertThrows(IllegalArgumentException.class, () -> vd.insertPointSites(bad));
		assertEquals(0, vd.numPointSites());

		var withNull = new ArrayList<Point>();
		withNull.add(new Point(1, 1));
		withNull.add(null);
		assertThrows(IllegalArgumentException.class, () -> vd.insertPointSites(withNull));
		assertEquals(0, vd.numPointSites());
		assertTrue(vd.check());
	}

	@Test
	public void newPointsRejectedAfterLineSitesButDuplicatesAllowed() {
		var vd = new VoronoiDiagram();
		var a = new Point(0, 0);
		var b = new Point(100, 0);
		var handles = vd.insertPointSites(List.of(a, b));
		vd.insertLineSite(handles.get(0), handles.get(1));

		// all-duplicates: allowed, returns existing handles, mutates nothing
		var again = vd.insertPointSites(List.of(a, b));
		assertSame(handles.get(0), again.get(0));
		assertSame(handles.get(1), again.get(1));
		assertEquals(2, vd.numPointSites());

		// any new point: rejected before mutation
		assertThrows(IllegalStateException.class, () -> vd.insertPointSites(List.of(a, new Point(50, 50))));
		assertEquals(2, vd.numPointSites());
	}

	@Test
	public void insertPolygonBuildsClosedBoundary() {
		var vd = new VoronoiDiagram();
		var square = List.of(new Point(-100, -100), new Point(100, -100), new Point(100, 100), new Point(-100, 100));
		var handles = vd.insertPolygon(square);
		assertEquals(4, handles.size());
		assertEquals(4, vd.numPointSites());
		assertEquals(4, vd.numLineSites());
		assertTrue(vd.check());
	}

	@Test
	public void insertPolygonMatchesManualPattern() {
		// jittered circle, bulk polygon vs manual points-then-segments
		var rnd = new Random(42);
		var n = 64;
		var poly = new ArrayList<Point>(n);
		for (var i = 0; i < n; i++) {
			var angle = 2 * Math.PI * i / n;
			var r = 300 * (1.0 + 0.1 * (rnd.nextDouble() - 0.5));
			poly.add(new Point(r * Math.cos(angle), r * Math.sin(angle)));
		}

		var manual = new VoronoiDiagram();
		var handles = new ArrayList<Vertex>(n);
		for (Point p : poly) {
			handles.add(manual.insertPointSite(p));
		}
		for (var i = 0; i < n; i++) {
			manual.insertLineSite(handles.get(i), handles.get((i + 1) % n));
		}

		var bulk = new VoronoiDiagram();
		bulk.insertPolygon(poly);

		assertEquals(manual.numPointSites(), bulk.numPointSites());
		assertEquals(manual.numLineSites(), bulk.numLineSites());
		assertEquals(manual.numVertices(), bulk.numVertices());
		assertEquals(manual.numFaces(), bulk.numFaces());
		assertTrue(manual.check());
		assertTrue(bulk.check());
	}

	@Test
	public void insertPolygonRejectsDegenerateInput() {
		var vd = new VoronoiDiagram();
		assertThrows(IllegalArgumentException.class, () -> vd.insertPolygon(null));
		assertThrows(IllegalArgumentException.class,
				() -> vd.insertPolygon(List.of(new Point(0, 0), new Point(1, 1))));
	}
}
