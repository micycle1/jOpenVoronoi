package org.rogach.jopenvoronoi.util;

/**
 * Maps 2D points to their index along a Hilbert space-filling curve.
 * <p>
 * Used to spatially sort point sites before bulk insertion: consecutive
 * Hilbert indices are spatially close, so kd-tree point location and the
 * delete-tree region stay cache-warm from one insertion to the next.
 */
public final class HilbertCurve {

	/** grid is 2^ORDER x 2^ORDER cells */
	public static final int ORDER = 16;

	private static final int GRID = 1 << ORDER;

	private HilbertCurve() {
		throw new AssertionError("HilbertCurve is a utility class");
	}

	/**
	 * Hilbert-curve index of {@code (x, y)} on a {@code 2^16 x 2^16} grid laid
	 * over the square {@code [min, max]^2}. Coordinates outside the square are
	 * clamped to the nearest cell.
	 *
	 * @return index in {@code [0, 2^32)}
	 */
	public static long index(double x, double y, double min, double max) {
		var scale = GRID / (max - min);
		var cx = clampToGrid((x - min) * scale);
		var cy = clampToGrid((y - min) * scale);
		return xy2d(cx, cy);
	}

	private static int clampToGrid(double v) {
		if (v <= 0) {
			return 0;
		}
		if (v >= GRID - 1) {
			return GRID - 1;
		}
		return (int) v;
	}

	/** classic Hilbert xy-to-distance conversion (Wikipedia's iterative form) */
	private static long xy2d(int x, int y) {
		var d = 0L;
		for (var s = GRID / 2; s > 0; s /= 2) {
			var rx = (x & s) > 0 ? 1 : 0;
			var ry = (y & s) > 0 ? 1 : 0;
			d += (long) s * s * ((3 * rx) ^ ry);
			// rotate quadrant
			if (ry == 0) {
				if (rx == 1) {
					x = s - 1 - x;
					y = s - 1 - y;
				}
				var t = x;
				x = y;
				y = t;
			}
		}
		return d;
	}
}
