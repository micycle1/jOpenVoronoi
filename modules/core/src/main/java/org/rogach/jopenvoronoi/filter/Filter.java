package org.rogach.jopenvoronoi.filter;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;

/**
 * Base class for Voronoi diagram filters.
 * <p>
 * Filtering in jOpenVoronoi is a post-processing step on an already constructed
 * {@link HalfEdgeDiagram}. Call {@link org.rogach.jopenvoronoi.VoronoiDiagram#filter(Filter)}
 * with one or more filter instances after inserting all sites. For each edge in
 * the diagram, the filter's {@link #apply(Edge)} method is evaluated and any
 * edge that returns {@code false} has its {@code valid} flag cleared.
 * <p>
 * Filters are cumulative: applying a second filter can only invalidate more
 * edges, because {@link org.rogach.jopenvoronoi.VoronoiDiagram#filter(Filter)}
 * does not reset previously invalidated edges. Call
 * {@link org.rogach.jopenvoronoi.VoronoiDiagram#filterReset()} before starting a
 * new filtering pass from the unfiltered diagram.
 * <p>
 * Downstream consumers such as {@link org.rogach.jopenvoronoi.offset.Offset}
 * inspect the {@code valid} flags to decide which parts of the diagram to use.
 */
public abstract class Filter {
	/** Diagram currently being filtered. */
	protected HalfEdgeDiagram g;

	/**
	 * Supplies the diagram being filtered before {@link #apply(Edge)} is called.
	 *
	 * @param g half-edge diagram owned by the enclosing Voronoi diagram
	 */
	public void setGraph(HalfEdgeDiagram g) {
		this.g = g;
	}

	/**
	 * Decides whether the given edge should remain valid.
	 *
	 * @param e edge being visited during a filtering pass
	 * @return {@code true} to keep the edge valid, {@code false} to clear its
	 *         {@code valid} flag
	 */
	public abstract boolean apply(Edge e);
};
