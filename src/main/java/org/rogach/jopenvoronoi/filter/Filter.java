package org.rogach.jopenvoronoi.filter;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;

/**
 * base-class for voronoi-diagram filters concrete sub-classes of Filter provide a predicate for determining if the edge belongs to the filtered graph.
 */
public abstract class Filter {
	/** vd-graph */
	protected HalfEdgeDiagram g;
	// set graph

	public void set_graph(HalfEdgeDiagram g) {
		this.g = g;
	}

	// does this edge belong to the filtered graph?
	public abstract boolean apply(Edge e);
};
