package org.rogach.jopenvoronoi.filter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class MedialAxisFilterTest {

	private static VoronoiDiagram square() {
		VoronoiDiagram vd = new VoronoiDiagram();
		List<Point> polygon = List.of(new Point(-1.0, -1.0), new Point(1.0, -1.0), new Point(1.0, 1.0), new Point(-1.0, 1.0));
		List<Vertex> vertexHandles = new ArrayList<>();
		for (Point point : polygon) {
			vertexHandles.add(vd.insertPointSite(point));
		}
		for (int i = 0; i < vertexHandles.size(); i++) {
			vd.insertLineSite(vertexHandles.get(i), vertexHandles.get((i + 1) % vertexHandles.size()));
		}
		return vd;
	}

	private static Set<Edge> medialAxisEdges(VoronoiDiagram vd) {
		return vd.getDiagram().edges.stream()
				.filter(edge -> edge.valid && edge.type != EdgeType.LINESITE && edge.type != EdgeType.NULLEDGE && edge.type != EdgeType.OUTEDGE)
				.collect(Collectors.toCollection(() -> Collections.newSetFromMap(new java.util.IdentityHashMap<>())));
	}

	@Test
	public void resultIsIndependentOfEdgeIterationOrder() {
		VoronoiDiagram baseline = square();
		baseline.filter(new PolygonInteriorFilter(true));
		baseline.filter(new MedialAxisFilter());
		Set<Point> baselinePositions = medialAxisEdges(baseline).stream().map(e -> e.source.position).collect(Collectors.toSet());

		// same construction, but the underlying edge list is shuffled before
		// filtering; a filter that lazily mutates/reads e.valid mid-pass could
		// produce a different pruning decision depending on this order.
		VoronoiDiagram shuffled = square();
		List<Edge> edges = shuffled.getDiagram().edges;
		Collections.shuffle(edges, new Random(12345));
		for (int i = 0; i < edges.size(); i++) {
			edges.get(i).diagramIndex = i;
		}
		shuffled.filter(new PolygonInteriorFilter(true));
		shuffled.filter(new MedialAxisFilter());
		Set<Point> shuffledPositions = medialAxisEdges(shuffled).stream().map(e -> e.source.position).collect(Collectors.toSet());

		Assertions.assertEquals(baselinePositions, shuffledPositions, "medial-axis filtering result must not depend on g.edges iteration order");
	}

	@Test
	public void filterInstanceIsReusableAcrossDiagrams() {
		MedialAxisFilter filter = new MedialAxisFilter();

		VoronoiDiagram first = square();
		first.filter(new PolygonInteriorFilter(true));
		first.filter(filter);
		int firstCount = medialAxisEdges(first).size();

		VoronoiDiagram second = square();
		second.filter(new PolygonInteriorFilter(true));
		second.filter(filter);
		int secondCount = medialAxisEdges(second).size();

		Assertions.assertEquals(firstCount, secondCount,
				"reusing a MedialAxisFilter instance on a second diagram must not return stale decisions from the first");
		Assertions.assertTrue(secondCount > 0);
	}

	@Test
	public void filterInstanceIsReusableAcrossRepeatedPasses() {
		MedialAxisFilter filter = new MedialAxisFilter();
		VoronoiDiagram vd = square();

		vd.filter(new PolygonInteriorFilter(true));
		vd.filter(filter);
		Set<Point> firstPass = medialAxisEdges(vd).stream().map(e -> e.source.position).collect(Collectors.toSet());

		// re-running the same filter instance on the same (unchanged) diagram,
		// after a reset, must reproduce exactly the same result: a filter that
		// accumulates per-edge state forever (rather than recomputing it fresh in
		// setGraph) could otherwise drift or short-circuit on the second pass.
		vd.filterReset();
		vd.filter(new PolygonInteriorFilter(true));
		vd.filter(filter);
		Set<Point> secondPass = medialAxisEdges(vd).stream().map(e -> e.source.position).collect(Collectors.toSet());

		Assertions.assertEquals(firstPass, secondPass, "reusing a MedialAxisFilter instance across repeated passes on the same diagram must be stable");
	}
}
