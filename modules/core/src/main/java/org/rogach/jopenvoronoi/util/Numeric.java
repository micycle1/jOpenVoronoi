package org.rogach.jopenvoronoi.util;

import java.util.ArrayList;
import java.util.List;

/**
 * Holds general numerical functions that are not specific to voronoi-diagrams
 * and may be useful elsewhere too
 */
public final class Numeric {

	private Numeric() {
		throw new AssertionError("Numeric is a utility class");
	}

	public static final double DEFAULT_CHOP_EPSILON = 1e-10;

	public static final double UNIT_INTERVAL_EPSILON = 1e-7;

	public static final double STRICT_ZERO_EPSILON = 1e-14;

	public static final double DOUBLE_COMPARISON_EPSILON = 1e-15;

	public static final double DISTANCE_EPSILON = 1e-3;

	public static final double SOLUTION_EDGE_EPSILON = 9e-4;

	public static final double BRENT_SOLVER_ABSOLUTE_EPSILON = 1e-20;

	// solve quadratic eqn: a*x*x + b*x + c = 0
	// returns real roots (0, 1, or 2) as vector
	public static List<Double> quadraticRoots(double a, double b, double c) {
		List<Double> roots = new ArrayList<>();
		if ((a == 0) && (b == 0)) {
			return roots;
		}
		if (a == 0) {
			roots.add(-c / b);
			return roots;
		}
		if (b == 0) {
			var sqr = -c / a;
			if (sqr > 0) {
				roots.add(Math.sqrt(sqr));
				roots.add(-roots.get(0));
				return roots;
			} else if (sqr == 0) {
				roots.add(0d);
				return roots;
			} else {
				// std::cout << " quadratic_roots() b == 0. no roots.\n";
				return roots;
			}
		}
		var disc = chop(b * b - 4 * a * c); // discriminant, chop!
		if (disc > 0) {
			double q;
			if (b > 0) {
				q = (b + Math.sqrt(disc)) / -2;
			} else {
				q = (b - Math.sqrt(disc)) / -2;
			}
			roots.add(q / a);
			roots.add(c / q);
			return roots;
		} else if (disc == 0) {
			roots.add(-b / (2 * a));
			return roots;
		}
		// std::cout << " quadratic_roots() disc < 0. no roots. disc= " << disc << "\n";
		return roots;
	}

	public static double determinant(double a, double b, double c, double d, double e, double f, double g, double h,
			double i) {
		return a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e);
	}

	public static double sq(double a) {
		return a * a;
	}

	public static double chop(double val, double tol) {
		if (Math.abs(val) < tol) {
			return 0;
		} else {
			return val;
		}
	}

	public static double chop(double val) {
		if (Math.abs(val) < DEFAULT_CHOP_EPSILON) {
			return 0;
		} else {
			return val;
		}
	}

	public static double snapUnitInterval(double val) {
		if (Math.abs(val) < UNIT_INTERVAL_EPSILON) {
			return 0.0;
		} else if (Math.abs(val - 1.0) < UNIT_INTERVAL_EPSILON) {
			return 1.0;
		}
		return val;
	}

	public static boolean areClose(double a, double b, double relativeTolerance, double absoluteTolerance) {
		var delta = Math.abs(a - b);
		if (delta < absoluteTolerance) {
			return true;
		}
		return delta <= relativeTolerance * Math.max(Math.abs(a), Math.abs(b));
	}

	public static double diangle(double x, double y) {
		if (y >= 0) {
			return (x >= 0 ? y / (x + y) : 1 - x / (-x + y));
		} else {
			return (x < 0 ? 2 - y / (-x - y) : 3 + x / (x - y));
		}
	}

	public static double diangleX(double a) {
		return (a < 2 ? 1 - a : a - 3);
	}

	public static double diangleY(double a) {
		return (a < 3 ? ((a > 1) ? 2 - a : a) : a - 4);
	}

	public static Pair<Double, Double> diangleXy(double a) {
		var x = diangleX(a);
		var y = diangleY(a);
		var norm = Math.sqrt(x * x + y * y);
		return new Pair<Double, Double>(x / norm, y / norm);
	}

	// return true if a lies in [less,more]
	public static boolean diangleBracket(double less, double a, double more) {
		if (less == more) {
			return false;
		} else if (less <= more) { // normal case..
			if ((less <= a) && (a < more)) {
				return true;
			}
		} else if (((less <= a) && (a <= 4)) || ((0 <= a) && (a < more))) {
			return true;
		}

		return false;
	}

	// return average of input angles
	public static double diangleMid(double alfa1, double alfa2) {
		if (alfa1 <= alfa2) {
			return (alfa1 + alfa2) / 2;
		} else {
			var opposite_mid = alfa2 + (alfa1 - alfa2) / 2;
			var mid = opposite_mid + 2;
			if (mid > 4) {
				mid = mid - 4;
			}
			assert ((0 <= mid) && (mid <= 4)) : " (0<=mid) && (mid<=4) ";
			return mid;
		}
	}

}
