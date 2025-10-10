package jopenvoronoi;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class SegmentsTest {

	@Test
	public void polygonSixVerticesTest() {
		VoronoiDiagram vd = new VoronoiDiagram();

		// create vertices as in the C++ example
		List<Point> pts = new ArrayList<>();
		pts.add(new Point(0.0, 0.5));
		pts.add(new Point(0.5, 0.5));
		pts.add(new Point(0.7, 0.0));
		pts.add(new Point(0.4, -0.3));
		pts.add(new Point(-0.2, -0.2));
		pts.add(new Point(0.0, 0.0));

		// insert point sites and collect handles
		List<Vertex> vertexHandles = new ArrayList<>();
		for (Point p : pts) {
			Vertex vh = vd.insert_point_site(p);
			vertexHandles.add(vh);
		}

		// now insert line segments connecting them into a closed polygon
		for (int n = 0; n < vertexHandles.size(); n++) {
			int next = (n == vertexHandles.size() - 1) ? 0 : n + 1;
			vd.insert_line_site(vertexHandles.get(n), vertexHandles.get(next));
		}

		// checks
		Assertions.assertTrue(vd.check(), "VoronoiDiagram.check() should return true");
		Assertions.assertEquals(6, vd.num_point_sites(), "Should have 6 point sites");
		Assertions.assertEquals(6, vd.num_line_sites(), "Should have 6 line sites");

		// get diagram and iterate faces/edges to exercise API
		HalfEdgeDiagram diagram = vd.getDiagram();
		Assertions.assertNotNull(diagram, "HalfEdgeDiagram should not be null");
		vd.getFaces().forEach(face -> {
			Point pos = face.site.position();
			Assertions.assertNotNull(pos, "face.site.position() should not be null");
			diagram.face_edges(face).forEach(edge -> {
				Assertions.assertNotNull(edge.source);
				Assertions.assertNotNull(edge.target);
				Assertions.assertNotNull(edge.source.position);
				Assertions.assertNotNull(edge.target.position);
			});
		});
	}

	@Test
	public void randomNonIntersectingSegmentsTest() {
		final int nmax = 30; // keep reasonably small for unit tests
		final double far = 1.0;
		// generate nmax random non-intersecting segments deterministically
		Random rnd = new Random(42);
		List<Segment> segs = new ArrayList<>();
		while (segs.size() < nmax) {
			Segment s;
			do {
				s = randomSegment(far, rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble());
			} while (segmentIntersects(segs, s));
			segs.add(s);
		}

		VoronoiDiagram vd = new VoronoiDiagram();
		List<Vertex[]> segmentHandles = new ArrayList<>(nmax);

		long t0 = System.nanoTime();
		for (Segment s : segs) {
			Vertex a = vd.insert_point_site(s.a);
			Vertex b = vd.insert_point_site(s.b);
			segmentHandles.add(new Vertex[] { a, b });
		}
		long tPointsNanos = System.nanoTime() - t0;

		// insert line sites
		long t1 = System.nanoTime();
		for (Vertex[] ids : segmentHandles) {
			vd.insert_line_site(ids[0], ids[1]);
		}
		long tLinesNanos = System.nanoTime() - t1;

		// basic checks
		Assertions.assertTrue(vd.check(), "VoronoiDiagram.check() should return true");
		Assertions.assertEquals(nmax, segmentHandles.size(), "Should have created nmax segments");
		Assertions.assertEquals(2 * nmax, vd.num_point_sites(), "Expect 2*nmax point sites (each segment had two endpoints)");
		Assertions.assertEquals(nmax, vd.num_line_sites(), "Should have nmax line sites");

		// print timings (informational)
		double tPointsSec = tPointsNanos / 1e9;
		double tLinesSec = tLinesNanos / 1e9;
		System.out.printf("Inserted %d point-sites in %.6f s, %d line-sites in %.6f s%n", 2 * nmax, tPointsSec, nmax, tLinesSec);

		// probe a few faces / vertices to ensure API works
		HalfEdgeDiagram diagram = vd.getDiagram();
		Assertions.assertNotNull(diagram, "HalfEdgeDiagram should not be null");
		// iterate small sample
		int checkedFaces = 0;
		for (Face f : vd.getFaces()) {
			if (checkedFaces++ > 10)
				break;
			Point pos = f.site.position();
			Assertions.assertNotNull(pos);
		}
	}

	// Helper small container for a segment (two Points)
	static class Segment {
		final Point a;
		final Point b;
	
		Segment(Point a, Point b) {
			this.a = a;
			this.b = b;
		}
	}

	// Vector cross product (2D) of vectors p and q (interpreted as position
	// vectors)
	private double cross(Point p, Point q) {
		return p.x * q.y - p.y * q.x;
	}

	// p2 - p1
	private Point sub(Point p2, Point p1) {
		return new Point(p2.x - p1.x, p2.y - p1.y);
	}

	// return true if segment s1 intersects with segment s2 (port of the C++ logic)
	private boolean intersects(Segment s1, Segment s2) {
		Point p1 = s1.a;
		Point p2 = s1.b;
		Point q1 = s2.a;
		Point q2 = s2.b;
		Point r = sub(p2, p1);
		Point s = sub(q2, q1);
		Point qp = sub(q1, p1);
	
		double rxs = cross(r, s);
		if (rxs == 0.0) { // parallel
			if (cross(qp, r) == 0.0) {
				// collinear => treat as intersection (as in C++ example)
				return true;
			} else {
				return false;
			}
		}
		double t = cross(qp, s) / rxs;
		double u = cross(qp, r) / rxs;
		return (0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0);
	}

	// test if s intersects with any of the segments in segs
	private boolean segmentIntersects(List<Segment> segs, Segment s) {
		for (Segment seg : segs) {
			if (intersects(seg, s))
				return true;
		}
		return false;
	}

	// create a random segment in the box determined by far (port of C++
	// random_segment)
	private Segment randomSegment(double far, double r1, double r2, double r3, double r4) {
		double pradius = (1.0 / Math.sqrt(2.0)) * far;
		double x1 = -pradius + 2 * pradius * r1;
		double y1 = -pradius + 2 * pradius * r2;
		double x2 = -pradius + 2 * pradius * r3;
		double y2 = -pradius + 2 * pradius * r4;
		return new Segment(new Point(x1, y1), new Point(x2, y2));
	}
}