package org.rogach.jopenvoronoi.filter;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Face;

/**
 * Retains only the part of the Voronoi diagram that lies inside a polygonal
 * boundary.
 * <p>
 * This filter is typically applied immediately after the polygon has been
 * inserted into the diagram and before a follow-up step such as
 * {@link org.rogach.jopenvoronoi.offset.Offset} generation or
 * {@link MedialAxisFilter} extraction. It keeps line-site and null edges, then
 * marks the remaining edges valid only on faces whose boundary orientation
 * matches the requested polygon interior side.
 * <p>
 * The constructor argument is expressed in terms of the direction used when
 * inserting each line site:
 * <ul>
 * <li>pass {@code true} when the polygon interior lies to the left of the
 * inserted segment direction (the common case for a counter-clockwise outer
 * contour)</li>
 * <li>pass {@code false} when the polygon interior lies to the right of the
 * inserted segment direction (the opposite winding)</li>
 * </ul>
 * Holes must therefore be inserted with the opposite winding of the outer
 * contour.
 */
public class PolygonInteriorFilter extends Filter {

	private final boolean interiorOnLeft;

	/**
	 * Creates a polygon interior filter for a known boundary orientation.
	 *
	 * @param interiorOnLeft {@code true} if the polygon interior is on the left
	 *                       side of each inserted boundary segment direction;
	 *                       {@code false} if it is on the right side
	 */
	public PolygonInteriorFilter(boolean interiorOnLeft) {
		this.interiorOnLeft = interiorOnLeft;
	}

	/**
	 * Keeps edges that belong to faces on the requested polygon-interior side.
	 */
	@Override
	public boolean apply(Edge e) {

		if (e.type == EdgeType.LINESITE || e.type == EdgeType.NULLEDGE) {
			return true;
		}

		var f = e.face;
		var s = f.getSite();
		if (s.isLine() && faceMatchesInteriorSide(f)) {
			return true;
		} else if (s.isPoint()) {
			var adjacentLineSite = findAdjacentLinesite(f);
			if (adjacentLineSite != null) {
				var twin = adjacentLineSite.twin;
				var twinFace = twin.face;
				if (faceMatchesInteriorSide(twinFace)) {
					return true;
				}
			} else {
				return false;
			}
		}
		return false;
	}

	/**
	 * Finds a face edge whose twin belongs to a neighboring line-site face.
	 *
	 * @param f point-site face to inspect
	 * @return an adjacent edge leading to a line-site face, or {@code null} when
	 *         no such edge exists
	 */
	private Edge findAdjacentLinesite(Face f) {
		var current = f.getEdge();
		var start = current;

		do {
			var twin = current.twin;
			if (twin != null) {
				var twf = twin.face;
				if (twf.getSite().isLine()) {
					return current;
				}
			}
			current = current.next;
		} while (current != start);

		return null;
	}

	/**
	 * Checks whether the line-site face corresponds to the requested interior
	 * side.
	 *
	 * @param f line-site face to inspect
	 * @return {@code true} when the face orientation matches the constructor
	 *         argument
	 */
	private boolean faceMatchesInteriorSide(Face f) {
		var current = f.getEdge();
		var start = current;
		do {
			if ((current.insertedDirection ? interiorOnLeft && current.type == EdgeType.LINESITE
					: !interiorOnLeft && current.type == EdgeType.LINESITE)) {
				return true;
			}
			current = current.next;
		} while (current != start);
		return false;
	}
}
