package org.rogach.jopenvoronoi.site;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class SiteApiTest {

	@Test
	public void exposesOnlyCamelCaseTypeChecks() throws NoSuchMethodException {
		Assertions.assertNotNull(Site.class.getMethod("isLine"));
		Assertions.assertNotNull(Site.class.getMethod("isPoint"));
		Assertions.assertNotNull(Site.class.getMethod("isArc"));
		Assertions.assertThrows(NoSuchMethodException.class, () -> Site.class.getMethod("is_linear"));
		Assertions.assertThrows(NoSuchMethodException.class, () -> Site.class.getMethod("is_quadratic"));
	}
}
