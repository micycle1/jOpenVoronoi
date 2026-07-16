package org.rogach.jopenvoronoi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Insertion-throughput benchmark harness. Run directly with {@code java}, not
 * via surefire: surefire enables assertions, and the per-insertion
 * {@code assert checker.isValid()} full-diagram walk would dominate the
 * measurement.
 *
 * <pre>
 * java -cp modules/core/target/classes;modules/core/target/test-classes org.rogach.jopenvoronoi.Perf [mode]
 *   (no args)        run the full suite: points 50k/100k + polygon 2000
 *   points  N seed   single points scenario
 *   polygon N seed   single jittered-circle polygon scenario (N segments)
 *   stress  N        one polygon build (for -ea stress runs), no warmup
 * </pre>
 *
 * Each scenario prints a diagram fingerprint: vertex/edge/face counts plus a
 * hash folded over the sorted {@code Double.doubleToLongBits} of all vertex
 * coordinates. Sorting makes the hash independent of graph-store ordering
 * (swap-remove reorders the vertex list); doubleToLongBits catches 1-ulp
 * drift. Fingerprints must stay identical across optimization steps.
 */
public class Perf {

	private static final int WARMUP_RUNS = 2;
	private static final int MEASURED_RUNS = 3;

	public static void main(String[] args) {
		warnIfAssertionsEnabled();
		if (args.length == 0) {
			runPoints(50_000, 123456);
			runPoints(50_000, 42);
			runPoints(100_000, 123456);
			runPolygon(2_000, 42);
		} else if (args[0].equals("points")) {
			runPoints(Integer.parseInt(args[1]), Long.parseLong(args[2]));
		} else if (args[0].equals("polygon")) {
			runPolygon(Integer.parseInt(args[1]), Long.parseLong(args[2]));
		} else if (args[0].equals("bulk")) {
			// bulk <uniform|clustered|xsorted> N seed : single-insert vs insertPointSites
			runBulkCompare(args[1], Integer.parseInt(args[2]), Long.parseLong(args[3]));
		} else if (args[0].equals("pointsloop")) {
			// repeated builds for profiling/steady-state: pointsloop N seed iterations
			int n = Integer.parseInt(args[1]);
			long seed = Long.parseLong(args[2]);
			int iterations = Integer.parseInt(args[3]);
			List<Point> pts = generateRandomPoints(n, seed);
			double best = 0;
			VoronoiDiagram vd = null;
			for (int i = 0; i < iterations; i++) {
				long t0 = System.nanoTime();
				vd = buildPoints(pts);
				long t1 = System.nanoTime();
				best = Math.max(best, n / ((t1 - t0) / 1e9));
			}
			System.out.printf("pointsloop n=%d seed=%d x%d: best %.1f sites/s%n", n, seed, iterations, best);
			System.out.println("  " + fingerprint(vd));
		} else if (args[0].equals("polygonloop")) {
			// repeated builds for profiling: polygonloop N seed iterations
			int n = Integer.parseInt(args[1]);
			long seed = Long.parseLong(args[2]);
			int iterations = Integer.parseInt(args[3]);
			List<Point> poly = generateJitteredCircle(n, seed);
			PolygonResult r = null;
			double bestPoints = 0, bestLines = 0;
			for (int i = 0; i < iterations; i++) {
				r = buildPolygon(poly);
				bestPoints = Math.max(bestPoints, r.pointsPerSec);
				bestLines = Math.max(bestLines, r.linesPerSec);
			}
			System.out.printf("polygonloop n=%d seed=%d x%d: best points %.1f, lines %.1f sites/s%n", n, seed,
					iterations, bestPoints, bestLines);
			System.out.println("  " + fingerprint(r.vd));
		} else if (args[0].equals("stress")) {
			int n = Integer.parseInt(args[1]);
			System.out.printf("stress polygon n=%d seed=42%n", n);
			PolygonResult r = buildPolygon(generateJitteredCircle(n, 42));
			System.out.printf("  points: %.1f sites/s, lines: %.1f sites/s%n", r.pointsPerSec, r.linesPerSec);
			System.out.println("  " + fingerprint(r.vd));
		} else {
			throw new IllegalArgumentException("unknown mode: " + args[0]);
		}
	}

	private static void warnIfAssertionsEnabled() {
		boolean ea = false;
		assert ea = true;
		if (ea) {
			System.out.println("WARNING: assertions are ENABLED - throughput numbers will be meaningless");
		}
	}

	// ------------------------------------------------------------------
	// Scenario A: uniform random points
	// ------------------------------------------------------------------

	private static void runPoints(int n, long seed) {
		List<Point> pts = generateRandomPoints(n, seed);
		for (int i = 0; i < WARMUP_RUNS; i++) {
			buildPoints(pts);
		}
		double[] rates = new double[MEASURED_RUNS];
		VoronoiDiagram last = null;
		for (int i = 0; i < MEASURED_RUNS; i++) {
			long t0 = System.nanoTime();
			last = buildPoints(pts);
			long t1 = System.nanoTime();
			rates[i] = n / ((t1 - t0) / 1e9);
		}
		System.out.printf("points n=%d seed=%d: median %.1f sites/s (runs: %s)%n", n, seed, median(rates),
				formatRates(rates));
		System.out.println("  " + fingerprint(last));
	}

	private static VoronoiDiagram buildPoints(List<Point> pts) {
		VoronoiDiagram vd = new VoronoiDiagram();
		for (Point p : pts) {
			vd.insertPointSite(p);
		}
		return vd;
	}

	private static List<Point> generateRandomPoints(int n, long seed) {
		Random rnd = new Random(seed);
		final double range = 1000.0;
		List<Point> pts = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			double x = rnd.nextDouble() * range - (range / 2.0);
			double y = rnd.nextDouble() * range - (range / 2.0);
			pts.add(new Point(x, y));
		}
		return pts;
	}

	// ------------------------------------------------------------------
	// Bulk insertion: single-insert loop vs insertPointSites
	// ------------------------------------------------------------------

	private static void runBulkCompare(String kind, int n, long seed) {
		List<Point> pts;
		switch (kind) {
			case "uniform":
				pts = generateRandomPoints(n, seed);
				break;
			case "clustered":
				pts = generateClusteredPoints(n, seed);
				break;
			case "xsorted":
				pts = generateRandomPoints(n, seed);
				pts.sort((a, b) -> Double.compare(a.x, b.x));
				break;
			default:
				throw new IllegalArgumentException("unknown bulk kind: " + kind);
		}
		for (int i = 0; i < WARMUP_RUNS; i++) {
			buildPoints(pts);
			buildBulk(pts);
		}
		double[] singleRates = new double[MEASURED_RUNS];
		double[] bulkRates = new double[MEASURED_RUNS];
		VoronoiDiagram lastSingle = null, lastBulk = null;
		for (int i = 0; i < MEASURED_RUNS; i++) {
			long t0 = System.nanoTime();
			lastSingle = buildPoints(pts);
			long t1 = System.nanoTime();
			lastBulk = buildBulk(pts);
			long t2 = System.nanoTime();
			singleRates[i] = n / ((t1 - t0) / 1e9);
			bulkRates[i] = n / ((t2 - t1) / 1e9);
		}
		System.out.printf("bulk %s n=%d seed=%d: single %.1f sites/s, bulk %.1f sites/s (%.2fx)%n", kind, n, seed,
				median(singleRates), median(bulkRates), median(bulkRates) / median(singleRates));
		System.out.println("  single " + fingerprint(lastSingle));
		System.out.println("  bulk   " + fingerprint(lastBulk));
	}

	private static VoronoiDiagram buildBulk(List<Point> pts) {
		VoronoiDiagram vd = new VoronoiDiagram();
		vd.insertPointSites(pts);
		return vd;
	}

	private static List<Point> generateClusteredPoints(int n, long seed) {
		Random rnd = new Random(seed);
		final int clusters = 100;
		final double range = 1000.0;
		final double sigma = 2.0;
		double[] cx = new double[clusters];
		double[] cy = new double[clusters];
		for (int i = 0; i < clusters; i++) {
			cx[i] = rnd.nextDouble() * range - (range / 2.0);
			cy[i] = rnd.nextDouble() * range - (range / 2.0);
		}
		List<Point> pts = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			int c = rnd.nextInt(clusters);
			pts.add(new Point(cx[c] + sigma * rnd.nextGaussian(), cy[c] + sigma * rnd.nextGaussian()));
		}
		return pts;
	}

	// ------------------------------------------------------------------
	// Scenario B: jittered-circle polygon (point phase + line phase)
	// ------------------------------------------------------------------

	private static void runPolygon(int segments, long seed) {
		List<Point> poly = generateJitteredCircle(segments, seed);
		for (int i = 0; i < WARMUP_RUNS; i++) {
			buildPolygon(poly);
		}
		double[] pointRates = new double[MEASURED_RUNS];
		double[] lineRates = new double[MEASURED_RUNS];
		PolygonResult last = null;
		for (int i = 0; i < MEASURED_RUNS; i++) {
			last = buildPolygon(poly);
			pointRates[i] = last.pointsPerSec;
			lineRates[i] = last.linesPerSec;
		}
		System.out.printf("polygon n=%d seed=%d: median points %.1f sites/s, lines %.1f sites/s%n", segments, seed,
				median(pointRates), median(lineRates));
		System.out.printf("  point runs: %s%n  line runs:  %s%n", formatRates(pointRates), formatRates(lineRates));
		System.out.println("  " + fingerprint(last.vd));
	}

	private static PolygonResult buildPolygon(List<Point> poly) {
		int n = poly.size();
		VoronoiDiagram vd = new VoronoiDiagram();
		List<Vertex> handles = new ArrayList<>(n);
		long t0 = System.nanoTime();
		for (Point p : poly) {
			handles.add(vd.insertPointSite(p));
		}
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			vd.insertLineSite(handles.get(i), handles.get((i + 1) % n));
		}
		long t2 = System.nanoTime();
		PolygonResult r = new PolygonResult();
		r.vd = vd;
		r.pointsPerSec = n / ((t1 - t0) / 1e9);
		r.linesPerSec = n / ((t2 - t1) / 1e9);
		return r;
	}

	private static class PolygonResult {
		VoronoiDiagram vd;
		double pointsPerSec;
		double linesPerSec;
	}

	private static List<Point> generateJitteredCircle(int segments, long seed) {
		Random rnd = new Random(seed);
		final double radius = 1000.0;
		final double jitter = 0.1;
		List<Point> pts = new ArrayList<>(segments);
		for (int i = 0; i < segments; i++) {
			double angle = 2 * Math.PI * i / segments;
			double r = radius * (1.0 + jitter * (rnd.nextDouble() - 0.5));
			pts.add(new Point(r * Math.cos(angle), r * Math.sin(angle)));
		}
		return pts;
	}

	// ------------------------------------------------------------------
	// Fingerprint + stats
	// ------------------------------------------------------------------

	static String fingerprint(VoronoiDiagram vd) {
		List<Vertex> vertices = vd.getVertices();
		long[] bits = new long[vertices.size() * 2];
		int i = 0;
		for (Vertex v : vertices) {
			bits[i++] = Double.doubleToLongBits(v.position.x);
			bits[i++] = Double.doubleToLongBits(v.position.y);
		}
		Arrays.sort(bits);
		long h = 1125899906842597L;
		for (long b : bits) {
			h = h * 31 + b;
		}
		return String.format("fingerprint: vertices=%d edges=%d faces=%d hash=%016x", vertices.size(),
				vd.getEdges().size(), vd.numAllFaces(), h);
	}

	private static double median(double[] values) {
		double[] sorted = values.clone();
		Arrays.sort(sorted);
		return sorted[sorted.length / 2];
	}

	private static String formatRates(double[] rates) {
		StringBuilder sb = new StringBuilder();
		for (double r : rates) {
			if (sb.length() > 0) {
				sb.append(", ");
			}
			sb.append(String.format("%.1f", r));
		}
		return sb.toString();
	}
}
