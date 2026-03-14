package org.rogach.jopenvoronoi;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.util.Pair;
import org.rogach.jopenvoronoi.vertex.Vertex;

/**
 * bundled BGL properties, see:
 * http:www.boost.org/doc/libs/1_44_0/libs/graph/doc/bundles.html
 *
 * dcel notes from http:www.holmes3d.net/graphics/dcel/
 * <p>
 * vertex (boost::out_edges) -leaving pointer to HalfEdge that has this vertex
 * as origin if many HalfEdges have this vertex as origin, choose one
 * arbitrarily
 * <p>
 * HalfEdge - origin pointer to vertex (boost::source) - face to the left of
 * halfedge - twin pointer to HalfEdge (on the right of this edge) - next
 * pointer to HalfEdge this edge starts from {@code h->twin->origin} and ends at
 * next vertex in {@code h->face} traveling ccw around boundary (allows face
 * traverse, follow {@code h->next} until we arrive back at h)
 * <p>
 * Face - edge pointer to HalfEdge this edge has this Face object as face
 * half-edge can be any one on the boundary of face special "infinite face",
 * face on "outside" of boundary may or may not store edge pointer
 * <p>
 * \brief half-edge diagram, based on the boost graph-library /
 * half_edge_diagram is a half-edge diagram class. Templated on Vertex/Edge/Face
 * property classes which allow attaching information to vertices/edges/faces
 * that is required for a particular algorithm.
 * <p>
 * Inherits from boost::adjacency_list minor additions allow storing
 * face-properties.
 * <p>
 * the hedi namespace contains functions for manipulating HEDIGraphs
 * <p>
 * For a general description of the half-edge data structure see e.g.: -
 * http:www.holmes3d.net/graphics/dcel/ - http:openmesh.org/index.php?id=228
 */
public class HalfEdgeDiagram {

	public Set<Vertex> vertices = new HashSet<>();
	public Set<Edge> edges = new HashSet<>();
	public Set<Face> faces = new HashSet<>();

	Vertex addVertex() {
		var v = new Vertex();
		vertices.add(v);
		return v;
	}

	// add a vertex with given properties, return vertex descriptor
	Vertex addVertex(Vertex v) {
		vertices.add(v);
		return v;
	}

	// return number of faces in graph
	public int numFaces() {
		return faces.size();
	}

	// return number of vertices in graph
	public int numVertices() {
		return vertices.size();
	}

	// return number of edges in graph
	public int numEdges() {
		return edges.size();
	}

	// return number of edges on Face f
	int numEdges(Face f) {
		return faceEdges(f).size();
	}

	// add an edge between vertices v1-v2
	Edge add_edge(Vertex v1, Vertex v2) {
		var e = new Edge(v1, v2);
		v1.outEdges.add(e);
		v2.inEdges.add(e);
		edges.add(e);
		return e;
	}

	// return true if v1-v2 edge exists
	boolean has_edge(Vertex v1, Vertex v2) {
		for (Edge e : v1.outEdges) {
			if (e.target == v2) {
				return true;
			}
		}
		return false;
	}

	// return v1-v2 Edge
	Edge edge(Vertex v1, Vertex v2) {
		for (Edge e : v1.outEdges) {
			if (e.target == v2) {
				return e;
			}
		}
		throw new RuntimeException("Edge not found in graph!");
	}

	// clear given vertex. this removes all edges connecting to the vertex.
	void clearVertex(Vertex v) {
		for (Edge e : v.outEdges) {
			e.target.inEdges.remove(e);
			edges.remove(e);
		}
		v.outEdges.clear();
		for (Edge e : v.inEdges) {
			e.source.outEdges.remove(e);
			edges.remove(e);
		}
		v.inEdges.clear();
	}

	// remove given vertex. call clear_vertex() before this!
	void removeVertex(Vertex v) {
		vertices.remove(v);
	}

	// remove given edge
	void removeEdge(Edge e) {
		e.source.outEdges.remove(e);
		e.target.inEdges.remove(e);
		edges.remove(e);
	}

	// delete a vertex. clear and remove.
	void deleteVertex(Vertex v) {
		clearVertex(v);
		removeVertex(v);
	}

	/**
	 * Insert vertex v into the middle of Edge e
	 *
	 * @param v vertex to insert
	 * @param e edge to split
	 */
	void addVertexInEdge(Vertex v, Edge e) {
		// the vertex v is inserted into the middle of edge e
		// edge e and its twin are replaced by four new edges: e1,e2 and their twins
		// te2,te1
		// before: face
		// e
		// previous-> source ------> target -> next
		// tw_next<- tw_trg <----- tw_src <- tw_previous
		// twin
		// twin_face
		//
		// after: face
		// e1 e2
		// previous-> source -> v -> target -> next
		// tw_next<- tw_trg <- v <- tw_src <- tw_previous
		// te2 te1
		// twin_face
		//

		var e_twin = e.twin;
		assert (e_twin != null) : " e_twin != null ";
		var esource = e.source;
		var etarget = e.target;
		var face = e.face;
		var twin_face = e_twin.face;
		var previous = previousEdge(e);
		var twin_previous = previousEdge(e_twin);

		assert (previous.face == e.face) : " previous.face == e.face ";
		assert (twin_previous.face == e_twin.face) : " twin_previous.face == e_twin.face ";

		var e1 = add_edge(esource, v);
		var te2 = add_edge(v, esource);
		twinEdges(e1, te2);

		var e2 = add_edge(v, etarget);
		var te1 = add_edge(etarget, v);
		twinEdges(e2, te1);

		// next-pointers
		previous.next = e1;
		e1.next = e2;
		e2.next = e.next;

		twin_previous.next = te1;
		te1.next = te2;
		te2.next = e_twin.next;

		// this copies params, face, k, type
		e1.copyFrom(e);
		e2.copyFrom(e);
		te1.copyFrom(e_twin);
		te2.copyFrom(e_twin);

		// update the faces
		face.edge = e1;
		twin_face.edge = te1;

		// finally, remove the old edge
		removeEdge(e);
		removeEdge(e_twin);
	}

	/**
	 * Adds two edges: one from v1 to v2, and one from v2 to v1
	 *
	 * @param v1 source vertex of the first edge
	 * @param v2 source vertex of the twin edge
	 * @return the two created twin edges
	 */
	Pair<Edge, Edge> addTwinEdges(Vertex v1, Vertex v2) {
		var e1 = add_edge(v1, v2);
		var e2 = add_edge(v2, v1);
		e1.twin = e2;
		e2.twin = e1;

		Edge canonical = chooseCanonicalByCoords(e1, e2);

		e1.base = canonical;
		e2.base = canonical;
		return new Pair<Edge, Edge>(e1, e2);
	}

	// make e1 the twin of e2 (and vice versa)
	void twinEdges(Edge e1, Edge e2) {
		assert (e1.target == e2.source) : "e1.target == e2.source";
		assert (e1.source == e2.target) : "e1.source == e2.target";
		e1.twin = e2;
		e2.twin = e1;

		Edge canonical = chooseCanonicalByCoords(e1, e2);
		e1.base = canonical;
		e2.base = canonical;
	}

	private Edge chooseCanonicalByCoords(Edge e1, Edge e2) {
		// prefer half-edge with target.x > source.x
		Point a1 = (e1.source != null) ? e1.source.position : null;
		Point b1 = (e1.target != null) ? e1.target.position : null;

		if (a1 != null && b1 != null) {
			int cmp = Double.compare(b1.x, a1.x);
			if (cmp > 0)
				return e1;
			if (cmp < 0)
				return e2;
			// x equal -> compare y
			cmp = Double.compare(b1.y, a1.y);
			if (cmp > 0)
				return e1;
			if (cmp < 0)
				return e2;
		}
		// fallback deterministic tie-break
		return (System.identityHashCode(e1) <= System.identityHashCode(e2)) ? e1 : e2;
	}

	// add a face, with given properties
	public Face addFace() {
		var f = new Face();
		f.attachDiagram(this);
		faces.add(f);
		return f;
	}

	// return all vertices adjecent to given vertex
	 List<Vertex> adjacentVertices(Vertex v) {
		List<Vertex> adj = new ArrayList<>();
		for (Edge e : v.outEdges) {
			adj.add(e.target);
		}
		return adj;
	}

	// return all vertices of given face
	 public List<Vertex> faceVertices(Face face) {
		List<Vertex> verts = new ArrayList<>();
		var startedge = face.edge; // the edge where we start
		var start_target = startedge.target;
		verts.add(start_target);

		var current = startedge.next;
		var count = 0;
		do {
			var current_target = current.target;
			verts.add(current_target);
			assert (current.face == current.next.face) : "current.face == current.next.face";
			current = current.next;
			if (count > 30000) {
				throw new AssertionError("count < 30000");
			}
			count++;
		} while (current != startedge);
		return verts;
	}

	// return edges of face f as a vector
	// NOTE: it is faster to write a do-while loop in client code than to call this
	// function!
	 public List<Edge> faceEdges(Face f) {
		var start_edge = f.edge;
		var current_edge = start_edge;
		List<Edge> out = new ArrayList<>();
		do {
			out.add(current_edge);
			current_edge = current_edge.next;
		} while (current_edge != start_edge);
		return out;
	}

	// return the previous edge. traverses all edges in face until previous found.
	 public Edge previousEdge(Edge e) {
		var previous = e.next;
		while (previous.next != e) {
			previous = previous.next;
		}
		return previous;
	}

	// return adjacent faces to the given vertex
	List<Face> adjacentFaces(Vertex q) {
		Set<Face> face_set = new HashSet<Face>();
		for (Edge e : q.outEdges) {
			face_set.add(e.face);
		}
		return new ArrayList<Face>(face_set);
	}

	// remove given v1-v2 edge
	void removeEdge(Vertex v1, Vertex v2) {
		removeEdge(edge(v1, v2));
	}

	// remove given v1-v2 edge and its twin
	void removeTwinEdges(Vertex v1, Vertex v2) {
		assert (has_edge(v1, v2)) : " has_edge(v1,v2) ";
		assert (has_edge(v2, v1)) : " has_edge(v2,v1) ";
		removeEdge(edge(v1, v2));
		removeEdge(edge(v2, v1));
	}

	// remove a degree-two Vertex from the middle of an Edge
	// preserve edge-properties (next, face, k)
	void removeDeg2Vertex(Vertex v) {
		// face1 e[1]
		// v1_prev -> v1 -> SPLIT -> v2 -> v2_next
		// v1_next <- v1 <- SPLIT <- v2 <- v2_prev
		// e[0] face2
		//
		// is replaced with a single edge:
		// face1
		// v1_prev -> v1 ----------> v2 -> v2_next
		// v1_next <- v1 <---------- v2 <- v2_prev
		// face2

		var v_edges = v.outEdges;
		assert (v_edges.size() == 2) : " v_edges.size() == 2";
		assert (v_edges.get(0).source == v && v_edges.get(1).source == v) : " v_edges.get(0).source == v && v_edges.get(1).source == v ";

		var v1 = v_edges.get(0).target;
		var v2 = v_edges.get(1).target;
		var v1_next = v_edges.get(0).next;
		var v1_prev = previousEdge(v_edges.get(0).twin);
		var v2_next = v_edges.get(1).next;
		var v2_prev = previousEdge(v_edges.get(1).twin);
		var face1 = v_edges.get(1).face;
		var face2 = v_edges.get(0).face;

		var twinEdges = addTwinEdges(v1, v2);
		var new1 = twinEdges.getFirst();
		var new2 = twinEdges.getSecond();
		setNext(new1, v2_next);
		setNext(new2, v1_next);
		setNext(v2_prev, new2);
		setNext(v1_prev, new1);
		face1.edge = new1;
		face2.edge = new2;
		new1.copyFrom(v_edges.get(1));
		new2.copyFrom(v_edges.get(0));
		removeTwinEdges(v, v1);
		removeTwinEdges(v, v2);
		removeVertex(v);
	}

	// set next-pointer of e1 to e2
	void setNext(Edge e1, Edge e2) {
		assert (e1.target == e2.source) : " e1.target == e2.source ";
		e1.next = e2;
	}

	// form a face from the edge-list:
	// e1->e2->...->e1
	// for all edges, set edge.face=f, and edge.k=k
	void setNextCycle(List<Edge> list, Face f, double k) {
		f.edge = list.get(0);
		for (var q = 0; q < list.size(); q++) {
			var e = list.get(q);
			e.face = f;
			e.k = k;
			if (q == list.size() - 1) {
				setNext(e, list.get(0));
			} else {
				setNext(e, list.get(q + 1));
			}
		}
	}

	// set next-pointers for the given list (but don't close to form a cycle)
	// also set face and k properties for edge
	void setNextChain(List<Edge> list, Face f, double k) {
		f.edge = list.get(0);
		for (var q = 0; q < list.size(); q++) {
			var e = list.get(q);
			e.face = f;
			e.k = k;
			if (q != list.size() - 1) {
				setNext(e, list.get(q + 1));
			}
		}
	}

	// set next-pointers for the list
	void setNextChain(List<Edge> list) {
		for (var q = 0; q < list.size(); q++) {
			var e = list.get(q);
			if (q != list.size() - 1) {
				setNext(e, list.get(q + 1));
			}
		}
	}

	// on a face, search and return the left/right edge from endp
	Pair<Edge, Edge> findNextPrev(Face f, Vertex endp) {
		var current = f.edge;
		var start_edge = current;
		Edge next_edge = null;
		Edge prev_edge = null;
		do {
			var src = current.source;
			var trg = current.target;
			if (src == endp) {
				next_edge = current;
			}
			if (trg == endp) {
				prev_edge = current;
			}
			current = current.next;
		} while (current != start_edge);
		assert (next_edge != null) : " next_edge != null ";
		assert (prev_edge != null) : " prev_edge != null ";
		return new Pair<Edge, Edge>(next_edge, prev_edge);
	}
}
