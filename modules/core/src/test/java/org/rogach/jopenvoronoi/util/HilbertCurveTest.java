package org.rogach.jopenvoronoi.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;

import org.junit.jupiter.api.Test;

public class HilbertCurveTest {

	@Test
	public void indexIsWithinRange() {
		assertTrue(HilbertCurve.index(-5000, -5000, -5000, 5000) >= 0);
		assertTrue(HilbertCurve.index(5000, 5000, -5000, 5000) < (1L << (2 * HilbertCurve.ORDER)));
	}

	@Test
	public void outOfRangeCoordinatesAreClamped() {
		assertEquals(HilbertCurve.index(-99999, -99999, -5000, 5000), HilbertCurve.index(-5000, -5000, -5000, 5000));
		assertEquals(HilbertCurve.index(99999, 0, -5000, 5000), HilbertCurve.index(5000, 0, -5000, 5000));
	}

	@Test
	public void smallGridIsABijectiveTour() {
		// on an 8x8 sub-grid aligned with the Hilbert grid, all indices are distinct
		var seen = new HashSet<Long>();
		var cell = 10000.0 / (1 << HilbertCurve.ORDER); // one grid cell
		for (var i = 0; i < 8; i++) {
			for (var j = 0; j < 8; j++) {
				seen.add(HilbertCurve.index(-5000 + (i + 0.5) * cell, -5000 + (j + 0.5) * cell, -5000, 5000));
			}
		}
		assertEquals(64, seen.size());
	}

	@Test
	public void consecutiveIndicesAreSpatiallyAdjacent() {
		// walk the first 4096 indices of the curve on the full grid by brute-force
		// inverse: instead, check locality statistically — consecutive-index cells
		// on a coarse 64x64 grid must be neighbours (Hilbert property on any
		// aligned power-of-two coarsening)
		var order = 6;
		var n = 1 << order;
		var cellSpan = 10000.0 / n;
		// map each coarse cell to its Hilbert index and verify each index 0..n*n-1
		// is used once and consecutive indices are 4-neighbours
		var pos = new int[n * n][];
		for (var i = 0; i < n; i++) {
			for (var j = 0; j < n; j++) {
				var x = -5000 + (i + 0.5) * cellSpan;
				var y = -5000 + (j + 0.5) * cellSpan;
				var d = HilbertCurve.index(x, y, -5000, 5000);
				// coarsen: full-order index / (cells per coarse cell)
				var coarse = (int) (d / (1L << (2 * (HilbertCurve.ORDER - order))));
				assertTrue(coarse >= 0 && coarse < n * n);
				pos[coarse] = new int[] { i, j };
			}
		}
		for (var k = 1; k < n * n; k++) {
			var a = pos[k - 1];
			var b = pos[k];
			var manhattan = Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]);
			assertEquals(1, manhattan, "cells at consecutive Hilbert indices must be grid neighbours");
		}
	}
}
