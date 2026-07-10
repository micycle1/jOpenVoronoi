package ags.utils.dataStructures.trees;

import ags.utils.dataStructures.BinaryHeap;
import ags.utils.dataStructures.MaxHeap;
import ags.utils.dataStructures.MinHeap;

public class KdTree<T> extends KdNode<T> {

	public KdTree(int dimensions) {
		this(dimensions, 24);
	}

	public KdTree(int dimensions, int bucketCapacity) {
		super(dimensions, bucketCapacity);
	}

	public NearestNeighborIterator<T> getNearestNeighborIterator(double[] searchPoint, int maxPointsReturned,
			DistanceFunction distanceFunction) {
		return new NearestNeighborIterator<T>(this, searchPoint, maxPointsReturned, distanceFunction);
	}

	/**
	 * Finds the single nearest neighbor to the search point.
	 * <p>
	 * Specialized 1-NN search: unlike {@link #findNearestNeighbors} it needs no
	 * heap bookkeeping, using a depth-first descent with bounding-box pruning.
	 *
	 * @return the value of the nearest point, or {@code null} if the tree is empty
	 */
	public T findNearestNeighbor(double[] searchPoint, DistanceFunction distanceFunction) {
		if (size() == 0) {
			return null;
		}
		var search = new NearestSearch<T>(searchPoint, distanceFunction);
		search.search(this);
		return search.bestValue();
	}

	private static final class NearestSearch<T> {
		private final double[] searchPoint;
		private final DistanceFunction distanceFunction;
		private double bestDist = Double.POSITIVE_INFINITY;
		private Object bestValue;

		NearestSearch(double[] searchPoint, DistanceFunction distanceFunction) {
			this.searchPoint = searchPoint;
			this.distanceFunction = distanceFunction;
		}

		@SuppressWarnings("unchecked")
		T bestValue() {
			return (T) bestValue;
		}

		void search(KdNode<T> node) {
			if (node.isLeaf()) {
				for (var i = 0; i < node.size; i++) {
					var d = distanceFunction.distance(node.points[i], searchPoint);
					if (d < bestDist) {
						bestDist = d;
						bestValue = node.data[i];
					}
				}
				return;
			}
			KdNode<T> near, far;
			if (searchPoint[node.splitDimension] > node.splitValue) {
				near = node.right;
				far = node.left;
			} else {
				near = node.left;
				far = node.right;
			}
			if (near.size > 0) {
				search(near); // descend the nearer child first so bestDist tightens early
			}
			if (far.size > 0
					&& distanceFunction.distanceToRect(searchPoint, far.minBound, far.maxBound) < bestDist) {
				search(far);
			}
		}
	}

	public MaxHeap<T> findNearestNeighbors(double[] searchPoint, int maxPointsReturned,
			DistanceFunction distanceFunction) {
		var pendingPaths = new BinaryHeap.Min<KdNode<T>>();
		var evaluatedPoints = new BinaryHeap.Max<T>();
		var pointsRemaining = Math.min(maxPointsReturned, size());
		pendingPaths.offer(0, this);

		while (pendingPaths.size() > 0 && (evaluatedPoints.size() < pointsRemaining
				|| (pendingPaths.getMinKey() < evaluatedPoints.getMaxKey()))) {
			nearestNeighborSearchStep(pendingPaths, evaluatedPoints, pointsRemaining, distanceFunction, searchPoint);
		}

		return evaluatedPoints;
	}

	@SuppressWarnings("unchecked")
	protected static <T> void nearestNeighborSearchStep(MinHeap<KdNode<T>> pendingPaths, MaxHeap<T> evaluatedPoints,
			int desiredPoints, DistanceFunction distanceFunction, double[] searchPoint) {
		// If there are pending paths possibly closer than the nearest evaluated point,
		// check it out
		var cursor = pendingPaths.getMin();
		pendingPaths.removeMin();

		// Descend the tree, recording paths not taken
		while (!cursor.isLeaf()) {
			KdNode<T> pathNotTaken;
			if (searchPoint[cursor.splitDimension] > cursor.splitValue) {
				pathNotTaken = cursor.left;
				cursor = cursor.right;
			} else {
				pathNotTaken = cursor.right;
				cursor = cursor.left;
			}
			var otherDistance = distanceFunction.distanceToRect(searchPoint, pathNotTaken.minBound,
					pathNotTaken.maxBound);
			// Only add a path if we either need more points or it's closer than furthest
			// point on list so far
			if (evaluatedPoints.size() < desiredPoints || otherDistance <= evaluatedPoints.getMaxKey()) {
				pendingPaths.offer(otherDistance, pathNotTaken);
			}
		}

		if (cursor.singlePoint) {
			var nodeDistance = distanceFunction.distance(cursor.points[0], searchPoint);
			// Only add a point if either need more points or it's closer than furthest on
			// list so far
			if (evaluatedPoints.size() < desiredPoints || nodeDistance <= evaluatedPoints.getMaxKey()) {
				for (var i = 0; i < cursor.size(); i++) {
					var value = (T) cursor.data[i];

					// If we don't need any more, replace max
					if (evaluatedPoints.size() == desiredPoints) {
						evaluatedPoints.replaceMax(nodeDistance, value);
					} else {
						evaluatedPoints.offer(nodeDistance, value);
					}
				}
			}
		} else {
			// Add the points at the cursor
			for (var i = 0; i < cursor.size(); i++) {
				var point = cursor.points[i];
				var value = (T) cursor.data[i];
				var distance = distanceFunction.distance(point, searchPoint);
				// Only add a point if either need more points or it's closer than furthest on
				// list so far
				if (evaluatedPoints.size() < desiredPoints) {
					evaluatedPoints.offer(distance, value);
				} else if (distance < evaluatedPoints.getMaxKey()) {
					evaluatedPoints.replaceMax(distance, value);
				}
			}
		}
	}
}
