package org.rogach.jopenvoronoi.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Modifier;

import org.junit.jupiter.api.Test;

public class NumericTest {

	@Test
	public void numericIsUtilityClass() throws Exception {
		var constructor = Numeric.class.getDeclaredConstructor();
		assertTrue(Modifier.isPrivate(constructor.getModifiers()));
		constructor.setAccessible(true);

		var error = assertThrows(InvocationTargetException.class, constructor::newInstance);
		assertInstanceOf(AssertionError.class, error.getCause());
	}

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

	@Test
	public void areCloseSupportsAbsoluteAndRelativeTolerance() {
		assertTrue(Numeric.areClose(0.0, Numeric.DOUBLE_COMPARISON_EPSILON * 0.5, Numeric.DISTANCE_EPSILON,
				Numeric.DOUBLE_COMPARISON_EPSILON));
		assertTrue(Numeric.areClose(1.0, 1.0 + Numeric.DOUBLE_COMPARISON_EPSILON * 0.5,
				Numeric.DISTANCE_EPSILON, Numeric.DOUBLE_COMPARISON_EPSILON));
		assertTrue(Numeric.areClose(1000.0, 1001.0, Numeric.DISTANCE_EPSILON, Numeric.DOUBLE_COMPARISON_EPSILON));
		assertTrue(Numeric.areClose(1000.0, 1000.5, Numeric.DISTANCE_EPSILON, Numeric.DOUBLE_COMPARISON_EPSILON));
		assertFalse(Numeric.areClose(1000.0, 1002.0, Numeric.DISTANCE_EPSILON, Numeric.DOUBLE_COMPARISON_EPSILON));
	}
}
