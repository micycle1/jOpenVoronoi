package org.rogach.jopenvoronoi.vertex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.util.Numeric;

/**
 * A vertex in the voronoi diagram an object of this type is held in the
 * BGL-graph for each vertex.
 */
public class Vertex {

	/** global vertex count TODO hold this in hedigraph instead? */
	private static int count = 0;

	// A map of this type is used by VoronoiDiagramChecker to check that all
	// vertices
	// have the expected (correct) degree (i.e. number of edges)
	// map for checking topology correctness
	public static Map<VertexType, Integer> expectedDegree = new HashMap<>();
	static {
		expectedDegree.put(VertexType.OUTER, 4); // special outer vertices
		expectedDegree.put(VertexType.NORMAL, 6); // normal vertex in the graph
		expectedDegree.put(VertexType.POINTSITE, 0); // point site
		expectedDegree.put(VertexType.ENDPOINT, 6); // end-point of line or arc
		expectedDegree.put(VertexType.SEPPOINT, 6); // end-point of separator
		expectedDegree.put(VertexType.SPLIT, 4); // split point, to avoid loops in delete-tree
		expectedDegree.put(VertexType.APEX, 4); // apex point on quadratic bisector
	}

	public List<Edge> outEdges = new ArrayList<>();
	public List<Edge> inEdges = new ArrayList<>();

	/**
	 * vertex status. updated/changed during an incremental graph update
	 */
	public VertexStatus status;
	/**
	 * The type of the vertex. Never(?) changes
	 */
	public VertexType type;
	/** TODO what is this? remove? */
	public double maxError;
	/** flag for indicating whether vertex is in the vertexQueue */
	public boolean inQueue;
	/** the position of the vertex. */
	public Point position;
	/** the offset-direction {-1,+1} of this vertex to the newly inserted site. */
	public double k3;
	/** diangle for a null-vertex and separator handling. */
	public double alfa;
	/** if this is a null-face, a handle to the null-face */
	public Face nullFace;
	/** the face of this vertex, if the vertex is a point-site */
	public Face face;
	/** clearance-disk radius, i.e. the closest Site is at this distance */
	public double r;
	public int diagramIndex = -1;

	public Vertex() {
	}

	public int degree() {
		return outEdges.size() + inEdges.size();
	}

	// ctor with given status and type
	public Vertex(Point p, VertexStatus st, VertexType t) {
		init(p, st, t);
	}

	// ctor with initial apex Point
	public Vertex(Point p, VertexStatus st, VertexType t, Point initDist) {
		init(p, st, t, initDist);
	}

	// ctor with initial k3-value
	public Vertex(Point p, VertexStatus st, VertexType t, Point initDist, double lk3) {
		init(p, st, t, initDist, lk3);
	}

	// ctor with initial clearance-disk radius
	public Vertex(Point p, VertexStatus st, VertexType t, double init_radius) {
		init(p, st, t);
		r = init_radius;
	}

	// set index, increase count, initialize in_queue to false.
	private void init() {
		count++;
		inQueue = false;
		alfa = -1; // invalid/non-initialized alfa value
		nullFace = null;
		type = VertexType.NORMAL;
		face = null;
		maxError = 0;
	}

	// set position and status
	private void init(Point p, VertexStatus st) {
		init();
		position = p;
		status = st;
	}

	// set position, status and type
	private void init(Point p, VertexStatus st, VertexType t) {
		init(p, st);
		type = t;
	}

	// set position, status, type, and clearance-disk through givem apex-point
	private void init(Point p, VertexStatus st, VertexType t, Point initDist) {
		init(p, st, t);
		initDist(initDist);
	}

	// set position, status, type, clerance-disk radius, and k3-side
	private void init(Point p, VertexStatus st, VertexType t, Point initDist, double lk3) {
		init(p, st, t, initDist);
		k3 = lk3;
	}

	// set in_queue false, and status to ::UNDECIDED
	public void resetStatus() {
		inQueue = false;
		status = VertexStatus.UNDECIDED;
	}

	public void setAlfa(Point dir) {
		alfa = Numeric.diangle(dir.x, dir.y);
	}

	// initialize clerance-disk
	public void initDist(Point p) {
		r = dist(p);
	}

	// return distance to a point from this vertex
	public double dist(Point p) {
		return position.sub(p).norm();
	}

	// set clearance-disk to zero
	public void zeroDist() {
		r = 0;
	}

	// return clearance disk-radius
	public double dist() {
		return r;
	}

	// in-circle predicate
	public double inCircle(Point p) {
		return dist(p) - r;
	}

	// reset the index count
	public static void resetCount() {
		count = 0;
	}

	@Override
	public String toString() {
		return String.format("V(%s)", position);
	}
}
