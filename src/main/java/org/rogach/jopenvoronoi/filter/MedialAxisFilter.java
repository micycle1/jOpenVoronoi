package org.rogach.jopenvoronoi.filter;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.vertex.Vertex;
import org.rogach.jopenvoronoi.vertex.VertexType;

/**
 * Retains an approximate medial axis of a Voronoi diagram.
 * <p>
 * After another filter has limited the diagram to the region of interest, this
 * filter keeps edges that behave like medial-axis branches and clears the
 * remaining ones by returning {@code false} from {@link #apply(Edge)}.
 */
public class MedialAxisFilter extends Filter {

	/**
	 * A dot-product threshold in [0,1] for filtering out edges between nearly
	 * parallel LineSite segments
	 */
	double dotProductThreshold;

	public MedialAxisFilter() {
		dotProductThreshold = 0.8;
	}

	/**
	 * @param threshold dot-product threshold in {@code [0,1]} for filtering out
	 *                  edges between nearly parallel line-site segments
	 */
	public MedialAxisFilter(double threshold) {
		dotProductThreshold = threshold;
	}

	/**
	 * Keeps edges that belong to the approximate medial axis.
	 */
	@Override
	public boolean apply(Edge e) {
		if (e.type == EdgeType.LINESITE || e.type == EdgeType.NULLEDGE) {
			return true;
		}
		if (e.type == EdgeType.SEPARATOR) {
			return false;
		}

		if (bothEndpointsPositive(e)) {
			return true;
		}

		if (segmentsParallel(e)) {
			return false;
		}

		return true;
	}

	/**
	 * @return {@code true} when both edge endpoints have a positive clearance
	 *         radius
	 */
	private boolean bothEndpointsPositive(Edge e) {
		var src = e.source;
		var trg = e.target;
		return (src.dist() > 0) && (trg.dist() > 0);
	}

	/**
	 * @return {@code true} when the adjacent polygon segments are nearly parallel
	 */
	private boolean segmentsParallel(Edge e) {
		var endp1 = findEndpoint(e);
		var endp2 = findEndpoint(e.twin);
		var e1 = findSegment(endp1);
		var e2 = findSegment(endp2);
		e2 = e2.twin;
		var dotprod = edgeDotprod(e1, e2);
		return dotprod > dotProductThreshold;
	}

	/**
	 * Calculate the dot-product between unit vectors aligned along edges
	 * {@code e1 -> e2}.
	 * <p>
	 * Since {@code e1} and {@code e2} are both line-sites, the direction is
	 * determined directly from their endpoints. Support for arc-site tangents
	 * would require additional logic here.
	 */
	private double edgeDotprod(Edge e1, Edge e2) {
		var src1 = e1.source;
		var trg1 = e1.target;
		var src2 = e2.source;
		var trg2 = e2.target;
		var sp1 = src1.position;
		var tp1 = trg1.position;
		var sp2 = src2.position;
		var tp2 = trg2.position;

		var dir1 = tp1.sub(sp1);
		var dir2 = tp2.sub(sp2);
		dir1.normalize();
		dir2.normalize();
		return dir1.dot(dir2);
	}

	/**
	 * Finds the line-site edge incident to the given endpoint vertex.
	 *
	 * @param v endpoint vertex on the polygon boundary
	 * @return incident line-site edge
	 */
	Edge findSegment(Vertex v) {
		for (Edge e : v.outEdges) {
			if (e.type == EdgeType.LINESITE) {
				return e;
			}
		}
		throw new RuntimeException("Failed to find line segment from vertex");
	}

	/**
	 * Finds the endpoint vertex reached through a neighboring null edge.
	 *
	 * @param e edge from the Voronoi diagram
	 * @return polygon endpoint associated with {@code e}
	 */
	Vertex findEndpoint(Edge e) {
		var next = e.next;
		var prev = g.previousEdge(e);
		Vertex endp;
		if (next.type == EdgeType.NULLEDGE) {
			endp = next.target;
			assert (endp.type == VertexType.ENDPOINT) : "endp.type == VertexType.ENDPOINT ";
		} else if (prev.type == EdgeType.NULLEDGE) {
			endp = prev.source;
			assert (endp.type == VertexType.ENDPOINT) : "endp.type == VertexType.ENDPOINT ";
		} else {
			throw new RuntimeException("Failed to find endpoint");
		}
		return endp;
	}

};
