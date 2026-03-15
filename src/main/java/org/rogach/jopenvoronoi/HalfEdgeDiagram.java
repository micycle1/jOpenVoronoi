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

	public List<Vertex> vertices = new ArrayList<>();
	public List<Edge> edges = new ArrayList<>();
	public List<Face> faces = new ArrayList<>();

	Vertex addVertex() {
	    var v = new Vertex();
	    registerVertex(v);
	    return v;
	}

	Vertex addVertex(Vertex v) {
	    registerVertex(v);
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
	    registerEdge(e);
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
	    while (!v.outEdges.isEmpty()) {
	        removeEdge(v.outEdges.get(v.outEdges.size() - 1));
	    }
	    while (!v.inEdges.isEmpty()) {
	        removeEdge(v.inEdges.get(v.inEdges.size() - 1));
	    }
	}

	// remove given vertex. call clear_vertex() before this!
	void removeVertex(Vertex v) {
	    unregisterVertexFromStore(v);
	}

	// remove given edge
	void removeEdge(Edge e) {
	    if (e.diagramIndex < 0) {
	        throw new IllegalStateException("removeEdge on unregistered edge: " + e);
	    }
	    e.source.outEdges.remove(e);
	    e.target.inEdges.remove(e);
	    unregisterEdgeFromStore(e);
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
		// Split:
		// e : esource -> etarget
		// eTwin : etarget -> esource
		//
		// Into:
		// e : esource -> v (reused)
		// e2 : v -> etarget (new)
		// eTwin : etarget -> v (reused)
		// te2 : v -> esource (new)

		final Edge eTwin = e.twin;
		assert (eTwin != null) : "e.twin != null";

		final Vertex esource = e.source;
		final Vertex etarget = e.target;

		assert (eTwin.source == etarget) : "eTwin.source == etarget";
		assert (eTwin.target == esource) : "eTwin.target == esource";

		final Edge eNext = e.next;
		final Edge twinNext = eTwin.next;

		final Edge e2 = new Edge(v, etarget);
		final Edge te2 = new Edge(v, esource);

		copyEdgeData(e2, e);
		copyEdgeData(te2, eTwin);

		registerEdge(e2);
		registerEdge(te2);

		v.outEdges.add(e2);
		v.outEdges.add(te2);

		etarget.inEdges.add(e2);
		esource.inEdges.add(te2);

		etarget.inEdges.remove(e);
		e.target = v;
		v.inEdges.add(e);

		esource.inEdges.remove(eTwin);
		eTwin.target = v;
		v.inEdges.add(eTwin);

		setNext(e, e2);
		setNext(e2, eNext);

		setNext(eTwin, te2);
		setNext(te2, twinNext);

		e.twin = te2;
		te2.twin = e;

		e2.twin = eTwin;
		eTwin.twin = e2;

		Edge base1 = chooseCanonicalByCoords(e, te2);
		e.base = base1;
		te2.base = base1;

		Edge base2 = chooseCanonicalByCoords(e2, eTwin);
		e2.base = base2;
		eTwin.base = base2;
	}

	private static void copyEdgeData(Edge dst, Edge src) {
		dst.face = src.face;
		dst.nullFace = src.nullFace;
		dst.hasNullFace = src.hasNullFace;
		dst.k = src.k;
		dst.type = src.type;
		dst.valid = src.valid;
		dst.sign = src.sign;
		dst.insertedDirection = src.insertedDirection;

		System.arraycopy(src.x, 0, dst.x, 0, src.x.length);
		System.arraycopy(src.y, 0, dst.y, 0, src.y.length);
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
	    registerFace(f);
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

	public Edge previousEdge(Edge e) {
		assert (e.prev != null) : "e.prev != null";
		return e.prev;
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
		// Before:
		//
		// face1:
		// v1_prev -> (v1 -> v) -> (v -> v2) -> v2_next
		//
		// face2:
		// v2_prev -> (v2 -> v) -> (v -> v1) -> v1_next
		//
		// After:
		//
		// face1:
		// v1_prev -> (v1 -> v2) -> v2_next
		//
		// face2:
		// v2_prev -> (v2 -> v1) -> v1_next

		final var vEdges = v.outEdges;
		assert (vEdges.size() == 2) : "v.outEdges.size() == 2";

		final Edge e0 = vEdges.get(0); // v -> v1 on face2
		final Edge e1 = vEdges.get(1); // v -> v2 on face1

		final Edge t0 = e0.twin; // v1 -> v on face1
		final Edge t1 = e1.twin; // v2 -> v on face2

		assert (t0 != null) : "e0.twin != null";
		assert (t1 != null) : "e1.twin != null";

		final Vertex v1 = e0.target;
		final Vertex v2 = e1.target;

		final Face face1 = e1.face;
		final Face face2 = e0.face;

		final Edge v1Next = e0.next;
		final Edge v2Next = e1.next;

		final Edge v1Prev = previousEdge(t0);
		final Edge v2Prev = previousEdge(t1);

		v.inEdges.remove(t0);
		t0.target = v2;
		v2.inEdges.add(t0);

		v.inEdges.remove(t1);
		t1.target = v1;
		v1.inEdges.add(t1);

		setNext(v1Prev, t0);
		setNext(t0, v2Next);

		setNext(v2Prev, t1);
		setNext(t1, v1Next);

		copyEdgeData(t0, e1);
		copyEdgeData(t1, e0);

		t0.twin = t1;
		t1.twin = t0;
		Edge canonical = chooseCanonicalByCoords(t0, t1);
		t0.base = canonical;
		t1.base = canonical;

		face1.edge = t0;
		face2.edge = t1;

		v.outEdges.remove(e0);
		v.outEdges.remove(e1);

		v1.inEdges.remove(e0);
		unregisterEdgeFromStore(e0);

		v2.inEdges.remove(e1);
		unregisterEdgeFromStore(e1);

		assert v.outEdges.isEmpty() : "v.outEdges empty";
		assert v.inEdges.isEmpty() : "v.inEdges empty";

		removeVertex(v);
	}

	// set next-pointer of e1 to e2
	void setNext(Edge e1, Edge e2) {
		assert (e1.target == e2.source) : " e1.target == e2.source ";
		e1.next = e2;
		e2.prev = e1;
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
	
	private void registerEdge(Edge e) {
	    assert e.diagramIndex == -1 : "edge already registered";
	    e.diagramIndex = edges.size();
	    edges.add(e);
	}

	private void unregisterEdgeFromStore(Edge e) {
	    if (e.diagramIndex < 0) {
	        throw new IllegalStateException("Edge not registered in diagram: " + e);
	    }

	    int idx = e.diagramIndex;
	    int lastIdx = edges.size() - 1;
	    Edge last = edges.get(lastIdx);

	    if (idx != lastIdx) {
	        edges.set(idx, last);
	        last.diagramIndex = idx;
	    }
	    edges.remove(lastIdx);
	    e.diagramIndex = -1;
	}

	private void registerVertex(Vertex v) {
	    assert v.diagramIndex == -1 : "vertex already registered";
	    v.diagramIndex = vertices.size();
	    vertices.add(v);
	}

	private void unregisterVertexFromStore(Vertex v) {
	    if (v.diagramIndex < 0) {
	        throw new IllegalStateException("Vertex not registered in diagram: " + v);
	    }

	    int idx = v.diagramIndex;
	    int lastIdx = vertices.size() - 1;
	    Vertex last = vertices.get(lastIdx);

	    if (idx != lastIdx) {
	        vertices.set(idx, last);
	        last.diagramIndex = idx;
	    }
	    vertices.remove(lastIdx);
	    v.diagramIndex = -1;
	}

	private void registerFace(Face f) {
	    assert f.diagramIndex == -1 : "face already registered";
	    f.diagramIndex = faces.size();
	    faces.add(f);
	}
}
