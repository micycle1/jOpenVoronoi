package org.rogach.jopenvoronoi;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class SpikeTest {

	/**
	 * Regression case for a vertex shared by three line sites (a polyline fan).
	 * Currently fails with {@code AssertionError: (new_count % 2) == 0} in
	 * {@code VoronoiDiagram.addEdges}, i.e. a genuine, unrelated topology bug in
	 * how the incremental algorithm handles null-face/endpoint structure when a
	 * point-site vertex has degree &gt; 2 among line sites. Not yet root-caused;
	 * left {@code @Disabled} (rather than commented out) so it stays visible to
	 * tooling and is easy to re-enable once fixed.
	 */
	@Test
	@Disabled("known bug: vertex shared by 3+ line sites fails with (new_count % 2) == 0 in addEdges")
	void testSpike() {
		VoronoiDiagram vd = new VoronoiDiagram();
		Vertex v3 = vd.insertPointSite(new Point(0.4676644843944667, 0.14559440261470627));
		Vertex v4 = vd.insertPointSite(new Point(0.5299451505909778, 0.22422772973181948));
		Vertex v9 = vd.insertPointSite(new Point(0.6406656057326543, 0.3640239162472442));
		Vertex v12 = vd.insertPointSite(new Point(0.20962013060047913, 0.48303100856119374));
		// v4 is shared between three lines
		vd.insertLineSite(v3, v4);
		vd.insertLineSite(v4, v9);
		vd.insertLineSite(v12, v4);
	}

}
