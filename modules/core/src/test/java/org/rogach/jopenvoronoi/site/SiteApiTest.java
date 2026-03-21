package org.rogach.jopenvoronoi.site;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.rogach.jopenvoronoi.geometry.Point;

public class SiteApiTest {

	@Test
	public void exposesOnlyCamelCaseTypeChecks() throws NoSuchMethodException {
		Assertions.assertNotNull(Site.class.getMethod("isLine"));
		Assertions.assertNotNull(Site.class.getMethod("isPoint"));
		Assertions.assertNotNull(Site.class.getMethod("isArc"));
		Assertions.assertThrows(NoSuchMethodException.class, () -> Site.class.getMethod("is_linear"));
		Assertions.assertThrows(NoSuchMethodException.class, () -> Site.class.getMethod("is_quadratic"));
	}

	@Test
	public void lineSiteInRegionTSnapsNearUnitIntervalEndpoints() {
		var line = new LineSite(new Point(0, 0), new Point(10, 0), 1);

		assertEquals(0.0, line.inRegionT(new Point(1e-8, 0)), 0.0);
		assertEquals(1.0, line.inRegionT(new Point(10 - 1e-8, 0)), 0.0);
	}

	@Test
	public void arcSiteInRegionTSnapsNearUnitIntervalEndpoints() {
		var arc = new ArcSite(new Point(1, 0), new Point(0, 1), new Point(0, 0), false);

		assertEquals(0.0, arc.inRegionT(new Point(1, 1e-8)), 0.0);
		assertEquals(1.0, arc.inRegionT(new Point(1e-8, 1)), 0.0);
	}
}
