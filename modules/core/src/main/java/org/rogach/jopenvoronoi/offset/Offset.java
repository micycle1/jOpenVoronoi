package org.rogach.jopenvoronoi.offset;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;

/**
 * Generates offsets from a Voronoi diagram.
 * <p>
 * An offset is always a closed loop. The loop consists of offset elements from
 * each face that the loop visits. Each face is associated with a {@link org.rogach.jopenvoronoi.site.Site Site}, and the
 * offset element from:
 * <ul>
 * <li>a point site is a circular arc</li>
 * <li>a line site is a line</li>
 * <li>an arc site is a circular arc</li>
 * </ul>
 * This class produces offsets at the given offset distance on the entire
 * Voronoi diagram. To produce offsets only inside or outside a given geometry,
 * apply a {@link org.rogach.jopenvoronoi.filter.Filter filter} first. The filter
 * sets the {@code valid} property of edges so that offsets are not produced on
 * faces with one or more invalid edges.
 */
public class Offset {
	/** vd-graph */
	HalfEdgeDiagram g;
	Set<Face> remainingFaces = new HashSet<>();
	/** list of output offsets */
	List<OffsetLoop> offsetList;

	/**
	 * @param g Voronoi diagram graph
	 */
	public Offset(HalfEdgeDiagram g) {
		this.g = g;
	}

	// create offsets at offset distance \a t
	public List<OffsetLoop> offset(double t) {
		offsetList = new ArrayList<OffsetLoop>();
		setFlags(t);
		Face start;
		var c = 0;
		while ((start = findStartFace()) != null) { // while there are faces that still require offsets
			offsetLoopWalk(start, t); // start on the face, and do an offset loop
			if (c > 30000) {
				throw new AssertionError("c > 30000, hang in offset walk");
			}
			c++;
		}

		return offsetList;
	}

	// find a suitable start face
	private Face findStartFace() {
		if (!remainingFaces.isEmpty()) {
			return remainingFaces.iterator().next();
		} else {
			return null;
		}
	}

	// perform an offset walk at given distance \a t,
	// starting at the given face
	private void offsetLoopWalk(Face start, double t) {
		var out_in_mode = false;
		var start_edge = findNextOffsetEdge(start.getEdge(), t, out_in_mode); // the first edge on the start-face
		if (start_edge == null) {
			throw new IllegalStateException("No bracketing edge found on start face for t=" + t);
		}
		var current_edge = start_edge;
		var loop = new OffsetLoop(); // store the output in this loop
		loop.offsetDistance = t;
		loop.add(new OffsetVertex(current_edge.point(t), current_edge));
		do {
			out_in_mode = edgeMode(current_edge, t);
			// find the next edge
			var next_edge = findNextOffsetEdge(current_edge.next, t, out_in_mode);
			if (next_edge == null) {
				throw new IllegalStateException("Broken offset walk: no next bracketing edge for t=" + t);
			}
			var current_face = current_edge.face;
			loop.add(offsetElementFromFace(current_face, current_edge, next_edge, t));
			remainingFaces.remove(current_face); // although we may revisit current_face (if it is non-convex), it
													// seems safe to mark it "done" here.
			current_edge = next_edge.twin;
		} while (current_edge != start_edge);
		offsetList.add(loop); // append the created loop to the output
	}

	// return an offset-element corresponding to the current face
	private OffsetVertex offsetElementFromFace(Face current_face, Edge current_edge, Edge next_edge, double t) {
		var s = current_face.getSite();
		var o = s.offset(current_edge.point(t), next_edge.point(t)); // ask the Site for offset-geometry here.
		var cw = true;
		if (!s.isLine()) { // point and arc-sites produce arc-offsets, for which cw must be set.
			cw = findCw(o.start(), o.center(), o.end()); // figure out cw or ccw arcs?
		}
		// add offset to output
		return new OffsetVertex(next_edge.point(t), o.radius(), o.center(), cw, current_face, next_edge);
	}

	/**
	 * Determine which bracketing direction applies to this edge for the current
	 * offset distance.
	 */
	private boolean edgeMode(Edge e, double t) {
		var src = e.source;
		var trg = e.target;
		var src_r = src.dist();
		var trg_r = trg.dist();
		if ((src_r < t) && (t < trg_r)) {
			return true;
		} else if ((trg_r < t) && (t < src_r)) {
			return false;
		} else {
			assert (false) : "failed to determine edge mode";
			return false;
		}
	}

	// figure out cw or ccw for an arc
	private boolean findCw(Point start, Point center, Point end) {
		return center.is_right(start, end); // NOTE: this only works for arcs smaller than a half-circle !
	}

	/**
	 * Starting at {@code e}, find the next edge on the face that brackets
	 * {@code t}.
	 * <p>
	 * We can be in one of two modes. If {@code mode == false} then we are looking
	 * for an edge where {@code src_t < t < trg_t}. If {@code mode == true} we are
	 * looking for an edge where {@code trg_t < t < src_t}.
	 */
	private Edge findNextOffsetEdge(Edge e, double t, boolean mode) {
		var start = e;
		var current = start;
		do {
			var src = current.source;
			var trg = current.target;
			var src_r = src.dist();
			var trg_r = trg.dist();
			if (!mode && (src_r < t) && (t < trg_r)) {
				return current;
			} else if (mode && (trg_r < t) && (t < src_r)) {
				return current;
			}
			current = current.next;
		} while (current != start);
		return null;
	}

	// go through all faces and set flag=0 if the face requires an offset.
	private void setFlags(double t) {
		for (Face f : g.faces) {
			var start = f.getEdge();
			var current = start;
			do {
				var src = current.source;
				var trg = current.target;
				var src_r = src.dist();
				var trg_r = trg.dist();
				if (tBracket(src_r, trg_r, t)) {
					remainingFaces.add(f);
				}
				current = current.next;
			} while (current != start);
		}

		// again go through faces again, and set flag=1 if any edge on the face is
		// invalid
		// this is required because an upstream filter will set valid=false on some
		// edges,
		// but not all, on a face where we do not want offsets.
		for (Face f : g.faces) {
			var start = f.getEdge();
			var current = start;
			do {
				if (!current.valid) {
					remainingFaces.remove(f); // don't offset faces with invalid edges
				}
				current = current.next;
			} while (current != start);
		}
	}

	// is t in (a,b) ?
	private boolean tBracket(double a, double b, double t) {
		var min_t = Math.min(a, b);
		var max_t = Math.max(a, b);
		return ((min_t < t) && (t < max_t));
	}

}
