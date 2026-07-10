package org.rogach.jopenvoronoi.jts;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.Polygon;

public class JtsVoronoiTest {

	private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();

	private static Polygon rectangle() {
		return GEOMETRY_FACTORY.createPolygon(new Coordinate[] { new Coordinate(0, 0), new Coordinate(4, 0),
				new Coordinate(4, 1), new Coordinate(0, 1), new Coordinate(0, 0) });
	}

	private static Polygon lShape() {
		return GEOMETRY_FACTORY.createPolygon(new Coordinate[] { new Coordinate(0, 0), new Coordinate(4, 0),
				new Coordinate(4, 1), new Coordinate(1, 1), new Coordinate(1, 3), new Coordinate(0, 3),
				new Coordinate(0, 0) });
	}

	@Test
	public void rejectsNonPolygonalInput() {
		assertThrows(IllegalArgumentException.class, () -> new JtsVoronoi(GEOMETRY_FACTORY.createPoint(new Coordinate(0, 0))));
		assertThrows(IllegalArgumentException.class, () -> new JtsVoronoi(GEOMETRY_FACTORY.createPolygon()));
		assertThrows(IllegalArgumentException.class, () -> new JtsVoronoi(null));
	}

	@Test
	public void interiorCellsAreASubsetOfAllCells() {
		JtsVoronoi voronoi = new JtsVoronoi(rectangle());

		List<Polygon> all = voronoi.getCells();
		List<Polygon> interior = voronoi.getInteriorCells();

		assertFalse(interior.isEmpty());
		assertTrue(interior.size() < all.size(), "exterior cells must be excluded from the interior view");
	}

	@Test
	public void medialAxisOfRectangleStaysInsideAndSpansTheLongAxis() {
		Polygon rect = rectangle();
		JtsVoronoi voronoi = new JtsVoronoi(rect);

		MultiLineString axis = voronoi.getMedialAxis();
		assertFalse(axis.isEmpty());
		assertTrue(rect.buffer(1e-6).covers(axis), "medial axis must lie inside the polygon");

		// Full axis of a 4x1 rectangle: two corner "V"s joined by the center
		// segment from (0.5, 0.5) to (3.5, 0.5), total length 3 + 4*sqrt(0.5).
		assertEquals(3.0 + 4.0 * Math.sqrt(0.5), axis.getLength(), 1e-6);
	}

	@Test
	public void prunedMedialAxisIsNotLongerThanFullAxis() {
		JtsVoronoi voronoi = new JtsVoronoi(lShape());

		double full = voronoi.getMedialAxis().getLength();
		double pruned = voronoi.getMedialAxis(1.0).getLength();

		assertTrue(full > 0);
		assertTrue(pruned <= full);
	}

	@Test
	public void medialAxisCoverageTilesTheRectangle() {
		Polygon rect = rectangle();
		JtsVoronoi voronoi = new JtsVoronoi(rect);

		List<Polygon> coverage = voronoi.getMedialAxisCoverage();

		// The full axis (corner branches + center segment) cuts the rectangle
		// into two triangles and two trapezoids.
		assertEquals(4, coverage.size());

		double area = coverage.stream().mapToDouble(Polygon::getArea).sum();
		assertEquals(rect.getArea(), area, 1e-9);
		for (Polygon piece : coverage) {
			assertTrue(rect.buffer(1e-6).covers(piece), "coverage piece must lie inside the polygon");
		}
	}

	@Test
	public void medialAxisCoverageTilesTheLShape() {
		Polygon shape = lShape();
		JtsVoronoi voronoi = new JtsVoronoi(shape);

		List<Polygon> coverage = voronoi.getMedialAxisCoverage();

		assertTrue(coverage.size() > 1);
		double area = coverage.stream().mapToDouble(Polygon::getArea).sum();
		assertEquals(shape.getArea(), area, 1e-9);
	}

	@Test
	public void medialAxisQueriesAreRepeatable() {
		JtsVoronoi voronoi = new JtsVoronoi(lShape());

		double full1 = voronoi.getMedialAxis().getLength();
		double pruned = voronoi.getMedialAxis(1.0).getLength();
		double full2 = voronoi.getMedialAxis().getLength();

		assertEquals(full1, full2, 0.0, "filter state from a pruned query must not leak into later queries");
		assertTrue(pruned <= full1);
	}
}
