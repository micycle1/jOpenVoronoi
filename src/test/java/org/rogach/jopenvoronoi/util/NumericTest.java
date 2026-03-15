package org.rogach.jopenvoronoi.util;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

public class NumericTest {

	@Test
	public void chopUsesNamedDefaultEpsilon() {
		assertEquals(0.0, Numeric.chop(Numeric.DEFAULT_CHOP_EPSILON * 0.5), 0.0);
		assertEquals(Numeric.DEFAULT_CHOP_EPSILON * 2, Numeric.chop(Numeric.DEFAULT_CHOP_EPSILON * 2), 0.0);
	}

	@Test
	public void snapUnitIntervalUsesNamedTolerance() {
		assertEquals(0.0, Numeric.snapUnitInterval(Numeric.UNIT_INTERVAL_EPSILON * 0.5), 0.0);
		assertEquals(1.0, Numeric.snapUnitInterval(1.0 - Numeric.UNIT_INTERVAL_EPSILON * 0.5), 0.0);
		assertEquals(Numeric.UNIT_INTERVAL_EPSILON * 2, Numeric.snapUnitInterval(Numeric.UNIT_INTERVAL_EPSILON * 2), 0.0);
		assertEquals(1.0 - Numeric.UNIT_INTERVAL_EPSILON * 2,
				Numeric.snapUnitInterval(1.0 - Numeric.UNIT_INTERVAL_EPSILON * 2), 0.0);
	}
}
