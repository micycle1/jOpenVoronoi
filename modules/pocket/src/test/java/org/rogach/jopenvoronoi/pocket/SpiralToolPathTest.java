package org.rogach.jopenvoronoi.pocket;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.filter.MedialAxisFilter;
import org.rogach.jopenvoronoi.filter.PolygonInteriorFilter;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class SpiralToolPathTest {

	private static final File DEBUG_DIR = new File("target/spiral-debug");

	private static VoronoiDiagram buildPolygonDiagram(List<Point> polygon) {
		VoronoiDiagram vd = new VoronoiDiagram();
		List<Vertex> vs = new ArrayList<>();
		for (Point p : polygon) {
			vs.add(vd.insertPointSite(p));
		}
		for (int i = 0; i < vs.size(); i++) {
			vd.insertLineSite(vs.get(i), vs.get((i + 1) % vs.size()));
		}
		vd.filter(new PolygonInteriorFilter(true));
		vd.filter(new MedialAxisFilter());
		return vd;
	}

	private static List<Point> starPolygon(int spikes, double rOuter, double rInner) {
		List<Point> pts = new ArrayList<>();
		for (int i = 0; i < spikes * 2; i++) {
			double a = Math.PI * i / spikes;
			double r = (i % 2 == 0) ? rOuter : rInner;
			pts.add(new Point(r * Math.cos(a), r * Math.sin(a)));
		}
		return pts;
	}

	@Test
	public void concaveStarSpiralStaysInsideAndIsSmoother() throws IOException {
		List<Point> polygon = starPolygon(7, 0.7, 0.3);
		VoronoiDiagram vd = buildPolygonDiagram(polygon);

		SpiralToolPath sp = new SpiralToolPath(vd.getDiagram());
		sp.setStepOver(0.03);
		sp.run();

		List<List<Point>> smooth = sp.getToolPathComponents();
		List<List<Point>> raw = sp.getRawToolPathComponents();
		Assertions.assertFalse(smooth.isEmpty(), "should produce at least one component");
		Assertions.assertEquals(raw.size(), smooth.size());

		for (List<Point> component : smooth) {
			Assertions.assertTrue(component.size() > 10, "spiral should have many points");
			for (Point p : component) {
				Assertions.assertFalse(Double.isNaN(p.x) || Double.isNaN(p.y), "no NaN points");
				Assertions.assertTrue(insideOrOnBoundary(p, polygon, 1e-6), "point outside pocket: " + p + " (dist=" + distToBoundary(p, polygon) + ")");
			}
		}

		// smoothing must reduce total turning (jaggedness)
		double rawTurn = totalTurning(raw.get(0));
		double smoothTurn = totalTurning(smooth.get(0));
		Assertions.assertTrue(smoothTurn < rawTurn, "smoothed turning " + smoothTurn + " should be < raw " + rawTurn);

		// path must be continuous: no jumps larger than a few step-overs
		for (List<Point> component : smooth) {
			for (int i = 1; i < component.size(); i++) {
				double d = component.get(i - 1).distance(component.get(i));
				Assertions.assertTrue(d < 0.3, "unexpected jump of " + d + " at index " + i);
			}
		}

		render(polygon, raw.get(0), new File(DEBUG_DIR, "star-raw.png"));
		render(polygon, smooth.get(0), new File(DEBUG_DIR, "star-smooth.png"));
	}

	@Test
	public void concaveBlobIsFullyCovered() throws IOException {
		// leaf-like wobbly blob with pronounced concave notches
		List<Point> polygon = new ArrayList<>();
		int n = 240;
		for (int i = 0; i < n; i++) {
			double a = 2 * Math.PI * i / n;
			double r = 0.45 + 0.18 * Math.sin(3 * a) + 0.08 * Math.cos(7 * a);
			polygon.add(new Point(r * Math.cos(a), r * Math.sin(a)));
		}
		VoronoiDiagram vd = buildPolygonDiagram(polygon);

		double stepOver = 0.03;
		SpiralToolPath sp = new SpiralToolPath(vd.getDiagram());
		sp.setStepOver(stepOver);
		sp.run();

		Assertions.assertFalse(sp.getToolPathComponents().isEmpty());
		List<Point> path = new ArrayList<>();
		for (List<Point> component : sp.getToolPathComponents()) {
			for (Point p : component) {
				Assertions.assertTrue(insideOrOnBoundary(p, polygon, 1e-6), "point outside pocket: " + p);
			}
			path.addAll(component);
		}

		// coverage: every interior grid point must lie within ~stepOver of the path
		// (this fails if any concave lobe is skipped by the spiral)
		int grid = 40;
		for (int gy = 0; gy <= grid; gy++) {
			for (int gx = 0; gx <= grid; gx++) {
				Point q = new Point(-0.75 + 1.5 * gx / grid, -0.75 + 1.5 * gy / grid);
				if (!insideOrOnBoundary(q, polygon, -1e-9) || distToBoundary(q, polygon) < stepOver) {
					continue;
				}
				double best = Double.MAX_VALUE;
				for (Point p : path) {
					double d = q.distance(p);
					if (d < best) {
						best = d;
					}
				}
				Assertions.assertTrue(best <= stepOver * 1.2, "uncovered interior point " + q + " nearest path distance " + best);
			}
		}

		render(polygon, sp.getToolPathComponents().get(0), new File(DEBUG_DIR, "blob-smooth.png"));
		render(polygon, sp.getRawToolPathComponents().get(0), new File(DEBUG_DIR, "blob-raw.png"));
	}

	@Test
	public void unfilteredDiagramIsRejected() {
		List<Point> polygon = starPolygon(5, 0.6, 0.3);

		// no filters at all
		VoronoiDiagram vd1 = new VoronoiDiagram();
		List<Vertex> vs = new ArrayList<>();
		for (Point p : polygon) {
			vs.add(vd1.insertPointSite(p));
		}
		for (int i = 0; i < vs.size(); i++) {
			vd1.insertLineSite(vs.get(i), vs.get((i + 1) % vs.size()));
		}
		Assertions.assertThrows(IllegalStateException.class, () -> new SpiralToolPath(vd1.getDiagram()));

		// interior filter only (medial-axis filter missing)
		vd1.filter(new PolygonInteriorFilter(true));
		Assertions.assertThrows(IllegalStateException.class, () -> new SpiralToolPath(vd1.getDiagram()));

		// both filters: accepted
		vd1.filter(new MedialAxisFilter());
		Assertions.assertDoesNotThrow(() -> new SpiralToolPath(vd1.getDiagram()));
	}

	@Test
	public void smoothingDisabledReturnsRawPath() {
		List<Point> polygon = starPolygon(5, 0.6, 0.3);
		VoronoiDiagram vd = buildPolygonDiagram(polygon);

		SpiralToolPath sp = new SpiralToolPath(vd.getDiagram());
		sp.setStepOver(0.05);
		sp.setSmoothing(0.0);
		sp.run();

		Assertions.assertFalse(sp.getToolPathComponents().isEmpty());
		Assertions.assertEquals(sp.getRawToolPathComponents().get(0).size(), sp.getToolPathComponents().get(0).size());
		for (Point p : sp.getToolPathComponents().get(0)) {
			Assertions.assertTrue(insideOrOnBoundary(p, polygon, 1e-6), "raw point outside pocket: " + p);
		}
	}

	@Test
	public void squarePocketWorks() throws IOException {
		List<Point> polygon = List.of(new Point(-0.4, -0.4), new Point(0.4, -0.4), new Point(0.4, 0.4), new Point(-0.4, 0.4));
		VoronoiDiagram vd = buildPolygonDiagram(polygon);

		SpiralToolPath sp = new SpiralToolPath(vd.getDiagram());
		sp.setStepOver(0.05);
		sp.run();

		Assertions.assertFalse(sp.getToolPathComponents().isEmpty());
		for (Point p : sp.getToolPathComponents().get(0)) {
			Assertions.assertTrue(insideOrOnBoundary(p, polygon, 1e-6), "point outside pocket: " + p);
		}
		render(polygon, sp.getToolPathComponents().get(0), new File(DEBUG_DIR, "square-smooth.png"));
	}

	// -------------------------------------------------------------------------
	// helpers
	// -------------------------------------------------------------------------

	private static boolean insideOrOnBoundary(Point p, List<Point> polygon, double tol) {
		if (distToBoundary(p, polygon) <= tol) {
			return true;
		}
		// ray casting
		boolean inside = false;
		int n = polygon.size();
		for (int i = 0, j = n - 1; i < n; j = i++) {
			Point a = polygon.get(i);
			Point b = polygon.get(j);
			if ((a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x) {
				inside = !inside;
			}
		}
		return inside;
	}

	private static double distToBoundary(Point p, List<Point> polygon) {
		double best = Double.MAX_VALUE;
		int n = polygon.size();
		for (int i = 0; i < n; i++) {
			Point a = polygon.get(i);
			Point b = polygon.get((i + 1) % n);
			best = Math.min(best, distToSegment(p, a, b));
		}
		return best;
	}

	private static double distToSegment(Point p, Point a, Point b) {
		double dx = b.x - a.x;
		double dy = b.y - a.y;
		double len2 = dx * dx + dy * dy;
		double t = len2 <= 0 ? 0 : Math.max(0, Math.min(1, ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2));
		double px = a.x + t * dx - p.x;
		double py = a.y + t * dy - p.y;
		return Math.sqrt(px * px + py * py);
	}

	/** Sum of absolute direction changes along the polyline, in radians. */
	private static double totalTurning(List<Point> path) {
		double sum = 0;
		double prevAngle = Double.NaN;
		for (int i = 1; i < path.size(); i++) {
			Point a = path.get(i - 1);
			Point b = path.get(i);
			if (a.distance(b) < 1e-12) {
				continue;
			}
			double angle = Math.atan2(b.y - a.y, b.x - a.x);
			if (!Double.isNaN(prevAngle)) {
				double d = Math.abs(angle - prevAngle);
				if (d > Math.PI) {
					d = 2 * Math.PI - d;
				}
				sum += d;
			}
			prevAngle = angle;
		}
		return sum;
	}

	private static void render(List<Point> polygon, List<Point> path, File out) throws IOException {
		int size = 1200;
		double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE, maxX = -Double.MAX_VALUE, maxY = -Double.MAX_VALUE;
		for (Point p : polygon) {
			minX = Math.min(minX, p.x);
			minY = Math.min(minY, p.y);
			maxX = Math.max(maxX, p.x);
			maxY = Math.max(maxY, p.y);
		}
		double span = Math.max(maxX - minX, maxY - minY);
		double margin = 0.05 * span;
		double scale = size / (span + 2 * margin);
		double ox = minX - margin;
		double oy = minY - margin;

		BufferedImage img = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = img.createGraphics();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, size, size);

		Path2D poly = new Path2D.Double();
		for (int i = 0; i < polygon.size(); i++) {
			Point p = polygon.get(i);
			double px = (p.x - ox) * scale;
			double py = size - (p.y - oy) * scale;
			if (i == 0) {
				poly.moveTo(px, py);
			} else {
				poly.lineTo(px, py);
			}
		}
		poly.closePath();
		g.setColor(new Color(230, 140, 20));
		g.setStroke(new BasicStroke(3f));
		g.draw(poly);

		g.setColor(Color.BLACK);
		g.setStroke(new BasicStroke(1.2f));
		for (int i = 1; i < path.size(); i++) {
			Point a = path.get(i - 1);
			Point b = path.get(i);
			g.draw(new Line2D.Double((a.x - ox) * scale, size - (a.y - oy) * scale, (b.x - ox) * scale, size - (b.y - oy) * scale));
		}
		g.dispose();

		out.getParentFile().mkdirs();
		ImageIO.write(img, "png", out);
	}
}
