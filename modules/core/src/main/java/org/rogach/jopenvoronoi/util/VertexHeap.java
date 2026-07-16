package org.rogach.jopenvoronoi.util;

import java.util.Arrays;

import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * Max-heap of vertices ordered by descending {@link Vertex#queueAbs}
 * (precomputed {@code Math.abs(queueScore)}).
 * <p>
 * The sift algorithms are transliterated from {@link java.util.PriorityQueue}
 * so that poll order — including the order of equal-priority vertices — is
 * identical to a {@code PriorityQueue} with the same comparator, but without
 * the comparator interface dispatch. Ordering of ties affects delete-tree
 * growth on degenerate inputs, so it must not be changed.
 */
public final class VertexHeap {

	private Vertex[] queue = new Vertex[256];
	private int size;

	/** the comparator: large |in_circle| first */
	private static int cmp(Vertex a, Vertex b) {
		return -Double.compare(a.queueAbs, b.queueAbs);
	}

	public boolean isEmpty() {
		return size == 0;
	}

	public void add(Vertex e) {
		var i = size;
		if (i >= queue.length) {
			queue = Arrays.copyOf(queue, queue.length * 2);
		}
		size = i + 1;
		if (i == 0) {
			queue[0] = e;
		} else {
			siftUp(i, e);
		}
	}

	public Vertex poll() {
		if (size == 0) {
			return null;
		}
		var result = queue[0];
		var s = --size;
		var x = queue[s];
		queue[s] = null;
		if (s != 0) {
			siftDown(0, x);
		}
		return result;
	}

	public void clear() {
		for (var i = 0; i < size; i++) {
			queue[i] = null;
		}
		size = 0;
	}

	private void siftUp(int k, Vertex x) {
		while (k > 0) {
			var parent = (k - 1) >>> 1;
			var e = queue[parent];
			if (cmp(x, e) >= 0) {
				break;
			}
			queue[k] = e;
			k = parent;
		}
		queue[k] = x;
	}

	private void siftDown(int k, Vertex x) {
		var half = size >>> 1;
		while (k < half) {
			var child = (k << 1) + 1;
			var c = queue[child];
			var right = child + 1;
			if (right < size && cmp(c, queue[right]) > 0) {
				child = right;
				c = queue[child];
			}
			if (cmp(x, c) <= 0) {
				break;
			}
			queue[k] = c;
			k = child;
		}
		queue[k] = x;
	}
}
