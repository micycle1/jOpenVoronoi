package ags.utils.dataStructures.trees.thirdGenKD;

/**
 *
 */
public class SquareEuclideanDistanceFunction implements DistanceFunction {

	@Override
	public double distance(double[] p1, double[] p2) {
		var d = 0D;

		for (var i = 0; i < p1.length; i++) {
			var diff = (p1[i] - p2[i]);
			d += diff * diff;
		}

		return d;
	}

	@Override
	public double distanceToRect(double[] point, double[] min, double[] max) {
		var d = 0D;

		for (var i = 0; i < point.length; i++) {
			var diff = 0D;
			if (point[i] > max[i]) {
				diff = (point[i] - max[i]);
			} else if (point[i] < min[i]) {
				diff = (point[i] - min[i]);
			}
			d += diff * diff;
		}

		return d;
	}
}