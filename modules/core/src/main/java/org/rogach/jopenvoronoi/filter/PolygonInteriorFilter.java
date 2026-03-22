package org.rogach.jopenvoronoi.filter;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Face;

import java.util.IdentityHashMap;

/**
 * Keeps Voronoi edges that belong to the polygon interior for a known boundary
 * winding convention.
 * <p>
 * Interior classification is contextual, so it is computed lazily and cached
 * per face inside this filter:
 * <ul>
 * <li>line-site faces are interior when their defining inserted segment has the
 * requested interior side,</li>
 * <li>point-site faces are interior when they are adjacent to at least one
 * interior line-site face.</li>
 * </ul>
 */
public class PolygonInteriorFilter extends Filter {

	private final boolean interiorOnLeft;
	private final IdentityHashMap<Face, Boolean> cache = new IdentityHashMap<>();

	/**
	 * Creates a polygon interior filter for a known boundary orientation.
	 *
	 * @param interiorOnLeft {@code true} for counter-clockwise outer boundaries
	 *                       (and clockwise holes), i.e. when the polygon interior
	 *                       lies to the left of each inserted boundary segment;
	 *                       {@code false} for the opposite winding
	 */
	public PolygonInteriorFilter(boolean interiorOnLeft) {
		this.interiorOnLeft = interiorOnLeft;
	}

	/**
	 * Returns whether the edge should be kept.
	 * <p>
	 * Line-site and null edges are always kept; all other edges are kept only if
	 * their incident face is classified as interior.
	 */
	@Override
	public boolean apply(Edge e) {
		return e.type == EdgeType.LINESITE || e.type == EdgeType.NULLEDGE || (e.face != null && isInterior(e.face));
	}

	/**
	 * Returns whether the given face lies on the polygon-interior side under the
	 * configured winding convention.
	 */
	private boolean isInterior(Face f) {
		return cache.computeIfAbsent(f, this::computeInterior);
	}

	/**
	 * Computes the interior classification for a face.
	 * <ul>
	 * <li>Null faces are never interior.</li>
	 * <li>Line-site faces are classified directly from their defining segment
	 * edge.</li>
	 * <li>Point-site faces inherit interior status if adjacent to an interior
	 * line-site face.</li>
	 * </ul>
	 */
	private boolean computeInterior(Face f) {
		if (f.isNullFace() || f.getSite() == null)
			return false;

		if (f.getSite().isLine()) {
			return lineFaceMatchesInteriorSide(f);
		}

		if (f.getSite().isPoint()) {
			for (Face adj : f.getAdjacentFaces()) {
				if (adj != null && !adj.isNullFace() && adj.getSite() != null && adj.getSite().isLine()) {
					if (isInterior(adj))
						return true;
				}
			}
		}

		return false;
	}

	/**
	 * Returns whether the given line-site face corresponds to the requested
	 * polygon-interior side.
	 */
	private boolean lineFaceMatchesInteriorSide(Face f) {
		Edge e = f.getSite().edge();
		if (e == null)
			return false;
		return e.insertedDirection ? interiorOnLeft : !interiorOnLeft;
	}
}