package org.rogach.jopenvoronoi.filter;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;

/**
 * Base class for Voronoi diagram filters.
 * <p>
 * Concrete subclasses of {@link Filter} provide a predicate for determining if
 * an {@link Edge} belongs to the filtered graph.
 */
public abstract class Filter {
	/** vd-graph */
	protected HalfEdgeDiagram g;
	// set graph

	public void setGraph(HalfEdgeDiagram g) {
		this.g = g;
	}

	// does this edge belong to the filtered graph?
	public abstract boolean apply(Edge e);
};
