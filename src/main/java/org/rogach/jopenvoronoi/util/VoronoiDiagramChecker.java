package org.rogach.jopenvoronoi.util;

import java.util.List;

import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.VoronoiDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.vertex.Vertex;
import org.rogach.jopenvoronoi.vertex.VertexStatus;

/**
 * Provides sanity-checks for the VoronoiDiagram class
 */
public class VoronoiDiagramChecker {

	/** vd-graph */
	private final HalfEdgeDiagram g;

	public VoronoiDiagramChecker(HalfEdgeDiagram gi) {
		this.g = gi;
	}
	
	public VoronoiDiagramChecker(VoronoiDiagram vd) {
		this(vd.getDiagram());
	}

	// overall sanity-check for the diagram, calls other sanity-check functions
	public boolean isValid() {
		return allFacesOk() && vertexDegreeOk();
	}

	// check that the diagram is of degree three.
	// however ::SPLIT and ::APEX vertices are of degree 2.
	private boolean vertexDegreeOk() {
		for (Vertex v : g.vertices) {
			if (v.degree() != Vertex.expectedDegree.get(v.type)) {
				return false;
			}
		}
		return true;
	}

	// check that all vertices in the input vector have status ::IN
	public static boolean verticesAllIN(List<Vertex> q) {
		for (Vertex v : q) {
			if (v.status != VertexStatus.IN) {
				return false;
			}
		}
		return true;
	}

	// check that all faces are ok. calls face_ok()
	private boolean allFacesOk() {
		for (Face f : g.faces) {
			if (!checkFace(f)) {
				return false;
			}
		}
		return true;
	}

	// check that the face is ok
	public static boolean checkFace(Face f) {
		var current_edge = f.getEdge();

		var start_edge = current_edge;
		var k = current_edge.k;
		if (!((k == 1) || (k == -1))) {
			return false;
		}
		if (f.getSite() != null) { // guard against null-faces that dont have Site
			if (f.getSite().isPoint()) {
				if (!(k == 1)) {
					return false;
				}
			}
		}
		var n = 0;
		do {
			if (current_edge.k != k) { // all edges should have the same k-value
				return false;
			}
			if (!currFaceEqualsNext(current_edge)) {// all edges should have the same face
				return false;
			}

			if (!checkEdge(current_edge)) {
				return false;
			}

			current_edge = current_edge.next;
			n++;
			assert (n < 10000) : " n < 10000 ";
		} while (current_edge != start_edge);
		return true;
	}

	// check that current edge and next-edge are on the same face
	private static boolean currFaceEqualsNext(Edge e) {
		if (e.face != e.next.face) {
			return false;
		}
		return true;
	}

	// sanity-check for edge
	public static boolean checkEdge(Edge e) {
		var src = e.source;
		var trg = e.target;
		var twine = e.twin;

		if (twine == null) {
			return true;
		} else if (!(e == twine.twin)) {
			return false;
		}

		var tw_src = twine.source;
		var tw_trg = twine.target;
		return ((src == tw_trg) && (trg == tw_src));
	}

}
