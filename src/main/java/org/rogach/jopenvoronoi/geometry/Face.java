package org.rogach.jopenvoronoi.geometry;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.vertex.Vertex;

public class Face {

	private Edge edge;
	private Site site;
	private FaceStatus status;
	private boolean is_null_face;
	private HalfEdgeDiagram diagram;
	public int diagramIndex = -1;

	public Face() {
	}

	/** Returns the representative/starting edge of this face. */
	public Edge getEdge() {
		return edge;
	}

	/** Sets the representative/starting edge of this face. */
	public void setEdge(Edge edge) {
		this.edge = edge;
	}

	public List<Edge> getEdges() {
		return requireDiagram().faceEdges(this);
	}

	public List<Vertex> getVertices() {
		return requireDiagram().faceVertices(this);
	}

	public List<Face> getAdjacentFaces() {
		Set<Face> adjacentFaces = new LinkedHashSet<>();
		for (Edge boundaryEdge : getEdges()) {
			if (boundaryEdge.twin != null && boundaryEdge.twin.face != null && boundaryEdge.twin.face != this) {
				adjacentFaces.add(boundaryEdge.twin.face);
			}
		}
		return new ArrayList<>(adjacentFaces);
	}

	public Site getSite() {
		return site;
	}

	/** Sets the generator site for this face. */
	public void setSite(Site site) {
		this.site = site;
	}

	public FaceStatus getStatus() {
		return status;
	}

	/** Sets the status of this face. */
	public void setStatus(FaceStatus status) {
		this.status = status;
	}

	/**
	 * Returns whether this face is an internal helper face rather than a
	 * user-visible Voronoi cell.
	 * <p>
	 * Null faces are introduced by the line-site insertion algorithm to model
	 * endpoint/null-edge structure in the half-edge graph.
	 */
	public boolean isNullFace() {
		return is_null_face;
	}

	/** Sets whether this face is a null face. */
	public void setNullFace(boolean isNullFace) {
		this.is_null_face = isNullFace;
	}

	public boolean isPointSiteFace() {
		return site != null && site.isPoint();
	}

	public boolean isLineSiteFace() {
		return site != null && site.isLine();
	}

	public boolean isArcSiteFace() {
		return site != null && site.isArc();
	}

	public void attachDiagram(HalfEdgeDiagram diagram) {
		this.diagram = diagram;
	}

	private HalfEdgeDiagram requireDiagram() {
		if (diagram == null) {
			throw new IllegalStateException("Face is not attached to a diagram");
		}
		return diagram;
	}

	@Override
	public String toString() {
		var sb = new StringBuilder();
		sb.append("F(");
		var current = edge;
		var c = 0;
		do {
			if (current == null) {
				break;
			}
			sb.append(current.source.position);
			sb.append(">");
			current = current.next;
			c++;
		} while (current != edge && c < 100);
		if (c >= 100) {
			sb.append("...");
		}
		sb.append(")");
		return sb.toString();
	}
}
