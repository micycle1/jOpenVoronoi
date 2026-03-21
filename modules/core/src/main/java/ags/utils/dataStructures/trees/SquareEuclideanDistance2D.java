package ags.utils.dataStructures.trees;

public class SquareEuclideanDistance2D implements DistanceFunction {

    @Override
    public double distance(double[] p1, double[] p2) {
        var dx = p1[0] - p2[0];
        var dy = p1[1] - p2[1];
        return dx * dx + dy * dy;
    }

    @Override
    public double distanceToRect(double[] point, double[] min, double[] max) {
        var dx = 0D;
        if (point[0] > max[0]) {
            dx = point[0] - max[0];
        } else if (point[0] < min[0]) {
            dx = point[0] - min[0];
        }

        var dy = 0D;
        if (point[1] > max[1]) {
            dy = point[1] - max[1];
        } else if (point[1] < min[1]) {
            dy = point[1] - min[1];
        }

        return dx * dx + dy * dy;
    }
}