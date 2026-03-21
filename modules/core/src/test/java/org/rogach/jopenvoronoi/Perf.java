//package org.rogach.jopenvoronoi;
//
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Random;
//
//import org.junit.jupiter.api.Assertions;
//import org.rogach.jopenvoronoi.geometry.Point;
//import org.rogach.jopenvoronoi.vertex.Vertex;
//
//public class Perf {
//
//	public static void main(String[] args) {
//		measureInsertRandomVertices(10000);
//
//	}
//	
//	public static void measureInsertRandomVertices(int N) {
//		// populate N random vertices first
//		List<Point> pts = new ArrayList<>(N);
//		Random rnd = new Random(123456); // fixed seed for reproducibility
//		final double range = 1000.0;
//		for (int i = 0; i < N; i++) {
//		    double x = rnd.nextDouble() * range - (range / 2.0);
//		    double y = rnd.nextDouble() * range - (range / 2.0);
//		    pts.add(new Point(x, y));
//		}
//
//		// prepare diagram and handle list
//		VoronoiDiagram vd = new VoronoiDiagram();
//		List<Vertex> vertexHandles = new ArrayList<>(N);
//
//		// measure insertion of point sites
//		long startNs = System.nanoTime();
//		for (Point p : pts) {
//		    Vertex vh = vd.insertPointSite(p);
//		    vertexHandles.add(vh);
//		}
//		long endNs = System.nanoTime();
//
//		long elapsedNs = endNs - startNs;
//		double elapsedMs = elapsedNs / 1_000_000.0;
//		double avgNs = (N > 0) ? ((double) elapsedNs / N) : 0.0;
//		double insertsPerSec = (elapsedNs > 0) ? (N / (elapsedNs / 1_000_000_000.0)) : Double.POSITIVE_INFINITY;
//
//		System.out.printf("Inserted %d points: total = (%.3f ms), avg = %.3f ns, throughput = %.1f inserts/sec%n",
//		        N, elapsedMs, avgNs, insertsPerSec);
//
//		// sanity checks (optional)
//		Assertions.assertTrue(vd.check(), "VoronoiDiagram.check() should return true after insertions");
//		Assertions.assertEquals(N, vd.numPointSites(), "Should have N point sites after insertion");
//		}
//
//}
