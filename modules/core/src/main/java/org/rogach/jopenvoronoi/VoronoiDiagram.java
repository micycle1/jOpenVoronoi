package org.rogach.jopenvoronoi;

import static org.rogach.jopenvoronoi.util.VoronoiDiagramChecker.checkEdge;
import static org.rogach.jopenvoronoi.util.VoronoiDiagramChecker.checkFace;
import static org.rogach.jopenvoronoi.util.VoronoiDiagramChecker.verticesAllIN;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import org.apache.commons.math3.analysis.solvers.AllowedSolution;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.rogach.jopenvoronoi.filter.Filter;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Face;
import org.rogach.jopenvoronoi.geometry.FaceStatus;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.LineSite;
import org.rogach.jopenvoronoi.site.PointSite;
import org.rogach.jopenvoronoi.site.Site;
import org.rogach.jopenvoronoi.util.Numeric;
import org.rogach.jopenvoronoi.util.Pair;
import org.rogach.jopenvoronoi.util.SplitPointError;
import org.rogach.jopenvoronoi.util.VoronoiDiagramChecker;
import org.rogach.jopenvoronoi.vertex.Vertex;
import org.rogach.jopenvoronoi.vertex.VertexPositioner;
import org.rogach.jopenvoronoi.vertex.VertexStatus;
import org.rogach.jopenvoronoi.vertex.VertexType;

import ags.utils.dataStructures.trees.DistanceFunction;
import ags.utils.dataStructures.trees.KdTree;
import ags.utils.dataStructures.trees.SquareEuclideanDistance2D;

/**
 * Incremental Voronoi diagram for point and line sites.
 * <p>
 * A typical workflow is to create a diagram, insert all point sites first, and
 * then insert line sites between those point-site vertices when needed.
 * Point sites must be added before any line sites.
 * <p>
 * The diagram is stored as a half-edge structure. Most callers will use the
 * public query methods on this class; {@link #getDiagram()} exposes the backing
 * graph for lower-level access.
 */
public class VoronoiDiagram {
	
	private static final DistanceFunction DISTANCE_FUNCTION = new SquareEuclideanDistance2D();

	// HELPER-CLASSES
	/**
	 * Helper class used for sanity-checking the diagram.
	 */
	protected VoronoiDiagramChecker vd_checker;
	/**
	 * kd-tree for nearest neighbour search during point Site insertion
	 */
	protected KdTree<PointSite> kdTree;
	/**
	 * an algorithm for positioning vertices
	 */
	protected VertexPositioner vpos;

	// DATA

	/**
	 * priority_queue for vertex for processing sorted by decreasing fabs() of
	 * in_circle-predicate, so that the vertices whose IN/OUT status we are 'most
	 * certain' about are processed first queue of vertices to be processed
	 */
	protected PriorityQueue<Pair<Vertex, Double>> vertexQueue = new PriorityQueue<>(1, new AbsComparison());

	/**
	 * the half-edge diagram of the vd
	 */
	protected HalfEdgeDiagram g = new HalfEdgeDiagram();
	/**
	 * sites must fall within a circle with radius farRadius
	 */
	protected double farRadius;
	/**
	 * number of point sites
	 */
	protected int numPsites;
	/**
	 * number of line-segment sites
	 */
	protected int numLsites;
	/**
	 * number of arc-sites
	 */
	protected int numAsites;
	/**
	 * temporary variable for INCIDENT faces, will be reset to NONINCIDENT after a
	 * site has been inserted
	 */
	protected List<Face> incidentFaces = new ArrayList<>();
	/**
	 * temporary variable for in-vertices, out-vertices that need to be reset after
	 * a site has been inserted
	 */
	protected Set<Vertex> modifiedVertices = new HashSet<>();
	/**
	 * IN-vertices of the delete-tree, i.e. vertices to be deleted during the
	 * current incremental site insertion.
	 */
	protected List<Vertex> inVertices = new ArrayList<>();

	/**
	 * Creates an empty Voronoi diagram with a default far radius of {@code 5000}.
	 */
	public VoronoiDiagram() {
		this(5000);
	}

	/**
	 * Creates an empty VoronoiDiagram
	 *
	 * @param farRadius the radius of a circle within which all sites must be
	 *                  located.
	 */
	public VoronoiDiagram(double farRadius) {
		kdTree = new KdTree<PointSite>(2);
		vd_checker = new VoronoiDiagramChecker(g); // helper-class that checks topology/geometry
		vpos = new VertexPositioner(g); // helper-class that positions vertices
		this.farRadius = farRadius;
		initialize();
		numPsites = 3;
		numLsites = 0;
		numAsites = 0;
		resetVertexCount();
	}

	/**
	 * Inserts a point site into the Voronoi diagram.
	 *
	 * @param x x-coordinate of the point site
	 * @param y y-coordinate of the point site
	 * @return handle to the inserted point site
	 */
	public Vertex insertPointSite(double x, double y) {
		return insertPointSite(new Point(x, y));
	}

	/**
	 * Insert a PointSite into the voronoi diagram.
	 *
	 * <p>
	 * All PointSites must be inserted before any LineSites or ArcSites are
	 * inserted.
	 * <p>
	 * Duplicate PointSites (i.e. points with the same x,y coordinates) resolve to
	 * the already inserted point-site vertex. All PointSites must be inserted
	 * before any LineSites or ArcSites are inserted.
	 *
	 * <p>
	 * This is roughly "algorithm A" from the Sugihara-Iri 1994 paper, page 15/50 /
	 * -# find the face that is closest to the new site, see FaceGrid -# among the
	 * vertices on the closest face, find the seed vertex, see {@link #findSeedVertex} -#
	 * grow the tree of IN-vertices, see {@link #augmentVertexSet} -# add new
	 * voronoi-vertices on all IN-OUT edges so they become IN-NEW-OUT, see
	 * {@link #addVertices} -# add new face by splitting each INCIDENT face into two parts
	 * by inserting a NEW-NEW edge. see {@link #addEdges} -# repair the next-pointers of
	 * faces that have been modified. see {@link #repairFace} -# remove IN-IN edges and
	 * IN-NEW edges, see {@link #removeVertexSet} -# reset vertex/face status to be ready
	 * for next incremental operation, see {@link #resetStatus}
	 *
	 * @param p position of site; must not be {@code null}
	 * @return integer handle to the inserted point. use this integer when inserting
	 *         lines/arcs with {@link #insertLineSite(Vertex, Vertex)}
	 */
	public Vertex insertPointSite(Point p) {
		if (p == null) {
			throw new IllegalArgumentException("Point site cannot be null.");
		}

		if (p.norm() >= farRadius) {
			throw new IllegalArgumentException(
					String.format("Point site %s lies outside the configured far radius %.3f.", p, farRadius));
		}
		assert (p.norm() < farRadius) : "p.norm() < farRadius";

		var nearest = kdTree.findNearestNeighbors(new double[] { p.x, p.y }, 1, DISTANCE_FUNCTION).getMax();

		if (nearest.position().equals(p)) {
			return nearest.vertex();
		}

		if (numLsites > 0) {
			throw new IllegalStateException(
					"Cannot insert point sites after line sites. Insert all point sites before inserting any line sites.");
		}

		numPsites++;

		var new_vert = g.addVertex(new Vertex(p, VertexStatus.OUT, VertexType.POINTSITE));
		var new_site = new PointSite(p);
		new_site.setVertex(new_vert);
		var v_seed = findSeedVertex(nearest.face, new_site);
		markVertex(v_seed, new_site);
		augmentVertexSet(new_site); // grow the tree to maximum size
		addVertices(new_site); // insert new vertices on IN-OUT edges
		var newface = addFace(new_site);
		new_vert.face = newface; // Vertices that correspond to point-sites have their .face property set!
		for (Face f : incidentFaces) {
			addEdges(newface, f); // add NEW-NEW edges on all INCIDENT faces
		}
		repairFace(newface);
		removeVertexSet(); // remove all IN vertices and adjacent edges
		resetStatus(); // reset all vertices to UNDECIDED

		assert (checkFace(newface)) : "face check failed for newly added face";
		assert (vd_checker.isValid()) : "diagram validity check failed";

		return new_vert;
	}

	/**
	 * Inserts a LineSite into the diagram
	 * <p>
	 * All PointSites must be inserted before any LineSites are inserted. All
	 * LineSites should be inserted before any ArcSitess are inserted. It is an
	 * error to insert a LineSite that intersects an existing LineSite in the
	 * diagram!
	 * <p>
	 * Ignored (no-op) if start and end are equal.
	 * <p>
	 * The basic idea of the incremental diagram update is similar to that in
	 * insertPointSite(). The major differences are: - the handling of null-faces
	 * at the endpoints of the LineSite. - addition of SEPARATOR edges - addition of
	 * SPLIT vertices during {@link #augmentVertexSet} - removal of SPLIT vertices at the
	 * end
	 * <p>
	 * The steps of the algorithm are: -# based on \a idx1 and \a idx2, look up the
	 * corresponding vertex descriptors. It is an error if these are not found. -#
	 * find a seed-vertex -# grow the delete-tree of IN vertices. -# create or
	 * modify the null-faces at the startpoint and endpoint of the LineSite -#
	 * create and add LINESITE pseudo-edges -# add NEW vertices on all IN-OUT edges.
	 * -# add up to four SEPARATOR edges, where applicable -# add non-separator
	 * edges by calling {@link #addEdges} on all INCIDENT faces -# repair the
	 * next-pointers of faces that have been modified. see {@link #repairFace} -# remove
	 * IN-IN edges and IN-NEW edges, see {@link #removeVertexSet} -# remove SPLIT
	 * vertices -# reset vertex/face status to be ready for next incremental
	 * operation, see {@link #resetStatus}
	 * 
	 * @param start startpoint of line-segment
	 * @param end   endpoint of line-segment
	 * @return {@code true} when the line site is inserted successfully
	 */
	public boolean insertLineSite(Vertex start, Vertex end) {
		if (start.equals(end)) {
			return false;
		}

		numLsites++;

		// find the vertices corresponding to idx1 and idx2
		start.status = VertexStatus.OUT;
		end.status = VertexStatus.OUT;
		start.zeroDist();
		end.zeroDist();

		// create a point which is left of src->trg
		// determine k (offset-dir) for this point
		// then we know which site/face is the k==+1 and which is k==-1
		var src_se = start.position;
		var trg_se = end.position;
		var left = src_se.add(trg_se).mult(0.5).add(trg_se.sub(src_se).xyPerp()); // this is used below and in
																					// find_null_face()

		var pos_site = new LineSite(end.position, start.position, +1);
		var neg_site = new LineSite(start.position, end.position, -1);

		var seed_face = start.face; // assumes this point-site has a face!

		// on the face of start-point, find the seed vertex
		var v_seed = findSeedVertex(seed_face, pos_site);
		markVertex(v_seed, pos_site);

		augmentVertexSet(pos_site); // it should not matter if we use pos_site or neg_site here
		// todo(?) sanity checks:
		// check that end_face is INCIDENT?
		// check that tree (i.e. inVertices) includes end_face_seed ?

		// process the null-faces here
		// returns new seg_start/end vertices, new or existing null-faces, and separator
		// endpoints (if separators should be added)
		var dir2 = start.position.sub(end.position);
		var dir1 = end.position.sub(start.position);

		var null_face_start = findNullFace(start, end, left, dir1, pos_site);
		var seg_start = null_face_start.segEndpoint;
		var start_null_face = null_face_start.nullFace;
		var pos_sep_start = null_face_start.posSepVertex;
		var neg_sep_start = null_face_start.negSepVertex;
		var start_to_null = null_face_start.faceToNull;

		var null_face_end = findNullFace(end, start, left, dir2, pos_site);
		var seg_end = null_face_end.segEndpoint;
		var end_null_face = null_face_end.nullFace;
		var pos_sep_end = null_face_end.posSepVertex;
		var neg_sep_end = null_face_end.negSepVertex;
		var end_to_null = null_face_end.faceToNull;

		// now safe to set the zero-face edge
		// in the collinear case, set the edge for the face that "disappears" to a null
		// edge
		if (start_to_null != null) {
			start_to_null.setEdge(start_null_face.getEdge());
		}
		if (end_to_null != null) {
			end_to_null.setEdge(end_null_face.getEdge());
		}

		// create LINESITE pseudo edges and faces
		var twinEdges = g.addTwinEdges(seg_end, seg_start);
		var pos_edge = twinEdges.getFirst();
		var neg_edge = twinEdges.getSecond();

		pos_edge.insertedDirection = false;
		neg_edge.insertedDirection = true;
		pos_edge.type = EdgeType.LINESITE;
		neg_edge.type = EdgeType.LINESITE;
		pos_edge.k = +1;
		neg_edge.k = -1;
		var pos_face = addFace(pos_site); // this face to the left of start->end edge
		var neg_face = addFace(neg_site); // this face is to the left of end->start edge
		pos_face.setEdge(pos_edge);
		neg_face.setEdge(neg_edge);
		pos_edge.face = pos_face;
		neg_edge.face = neg_face;

		// associate sites with LINESITE edges
		pos_site.setEdge(pos_edge);
		neg_site.setEdge(neg_edge);

		addVertices(pos_site); // add NEW vertices on all IN-OUT edges.

		// add SEPARATORS
		// find SEPARATOR targets first
		var pos_start_target = findSeparatorTarget(start.face, pos_sep_start);
		var neg_start_target = findSeparatorTarget(start.face, neg_sep_start);

		// add positive separator edge at start
		addSeparator(start.face, start_null_face, pos_start_target, pos_sep_start, pos_face.getSite(), neg_face.getSite());

		// add negative separator edge at start
		addSeparator(start.face, start_null_face, neg_start_target, neg_sep_start, pos_face.getSite(), neg_face.getSite());
		start.face.setStatus(FaceStatus.NONINCIDENT); // face is now done.
		assert (checkFace(start.face)) : "face check failed for start endpoint face";

		var pos_end_target = findSeparatorTarget(end.face, pos_sep_end);
		var neg_end_target = findSeparatorTarget(end.face, neg_sep_end);
		// add positive separator edge at end
		addSeparator(end.face, end_null_face, pos_end_target, pos_sep_end, pos_face.getSite(), neg_face.getSite());

		// add negative separator edge at end
		addSeparator(end.face, end_null_face, neg_end_target, neg_sep_end, pos_face.getSite(), neg_face.getSite());
		end.face.setStatus(FaceStatus.NONINCIDENT);
		assert (checkFace(end.face)) : "face check failed for end endpoint face";

		// add non-separator edges by calling addEdges on all INCIDENT faces
		for (Face f : incidentFaces) {
			if (f.getStatus() == FaceStatus.INCIDENT) {// end-point faces already dealt with in addSeparator()
				addEdges(pos_face, f, neg_face, new Pair<Vertex, Vertex>(seg_start, seg_end)); // each INCIDENT face is
																								// split into two parts:
																								// newface and f
			}
		}

		// new vertices and edges inserted. remove the delete-set, repair faces.

		removeVertexSet();

		repairFace(pos_face, new Pair<Vertex, Vertex>(seg_start, seg_end), new Pair<Face, Face>(start_to_null, end_to_null),
				new Pair<Face, Face>(start_null_face, end_null_face));
		assert (checkFace(pos_face)) : "face check failed for positive-side face";

		repairFace(neg_face, new Pair<Vertex, Vertex>(seg_start, seg_end), new Pair<Face, Face>(start_to_null, end_to_null),
				new Pair<Face, Face>(start_null_face, end_null_face));
		assert (checkFace(neg_face)) : "face check failed for negative-side face";

		// we are done and can remove split-vertices
		for (Face f : incidentFaces) {
			removeSplitVertex(f);
		}
		resetStatus();

		assert (checkFace(start_null_face)) : "face check failed for start null-face";
		assert (checkFace(end_null_face)) : "face check failed for end null-face";
		assert (checkFace(pos_face)) : "face check failed for positive-side face";
		assert (checkFace(neg_face)) : "face check failed for negative-side face";
		assert (vd_checker.isValid()) : "diagram validity check failed";

		return true;
	}

	/**
	 * Returns the configured far radius that bounds the valid site insertion area.
	 *
	 * @return the far radius used when initializing the diagram
	 */
	public double getFarRadius() {
		return farRadius;
	}

	/**
	 * Returns the user-visible Voronoi faces.
	 * <p>
	 * Internal null faces used as implementation scaffolding for line-site
	 * insertion are hidden from this view.
	 *
	 * @return all non-null faces in the diagram
	 */
	public List<Face> getFaces() {
		List<Face> faces = new ArrayList<>();
		for (Face face : g.faces) {
			if (!face.isNullFace()) {
				faces.add(face);
			}
		}
		return faces;
	}

	/**
	 * Returns every face in the underlying half-edge graph, including internal null
	 * faces used during line-site insertion.
	 *
	 * @return all faces in the underlying graph
	 */
	public List<Face> getAllFaces() {
		return new ArrayList<>(g.faces);
	}

	/**
	 * Returns the same user-visible faces as {@link #getFaces()}.
	 *
	 * @return all non-null faces in the diagram
	 */
	public List<Face> getNonNullFaces() {
		return getFaces();
	}

	/**
	 * Returns the face directly associated with a vertex.
	 * <p>
	 * Point-site vertices reference their ordinary face via {@link Vertex#face},
	 * while vertices on null-face structures reference {@link Vertex#nullFace}.
	 *
	 * @param vertex the vertex whose associated face should be returned
	 * @return the face associated with the vertex
	 * @throws IllegalArgumentException if {@code vertex} is {@code null} or is not
	 *                                  directly associated with a face
	 */
	public Face getFace(Vertex vertex) {
		if (vertex == null) {
			throw new IllegalArgumentException("vertex cannot be null");
		}
		if (vertex.face != null) {
			return vertex.face;
		}
		if (vertex.nullFace != null) {
			return vertex.nullFace;
		}
		throw new IllegalArgumentException("vertex is not directly associated with a face");
	}

	/**
	 * Returns the half-edges that bound a face.
	 *
	 * @param face the face whose boundary edges should be returned
	 * @return the boundary edges of {@code face}
	 * @throws IllegalArgumentException if {@code face} is {@code null}
	 */
	public List<Edge> getFaceEdges(Face face) {
		if (face == null) {
			throw new IllegalArgumentException("face cannot be null");
		}
		return face.getEdges();
	}

	/**
	 * Returns the vertices that lie on the boundary of a face.
	 *
	 * @param face the face whose boundary vertices should be returned
	 * @return the boundary vertices of {@code face}
	 * @throws IllegalArgumentException if {@code face} is {@code null}
	 */
	public List<Vertex> getFaceVertices(Face face) {
		if (face == null) {
			throw new IllegalArgumentException("face cannot be null");
		}
		return face.getVertices();
	}

	/**
	 * Returns the faces adjacent to a face across its boundary edges.
	 *
	 * @param face the face whose neighboring faces should be returned
	 * @return the faces adjacent to {@code face}
	 * @throws IllegalArgumentException if {@code face} is {@code null}
	 */
	public List<Face> getAdjacentFaces(Face face) {
		if (face == null) {
			throw new IllegalArgumentException("face cannot be null");
		}
		return face.getAdjacentFaces();
	}

	/**
	 * Returns all half-edges in the underlying diagram graph.
	 *
	 * @return every edge currently stored in the diagram
	 */
	public List<Edge> getEdges() {
		return new ArrayList<>(g.edges);
	}

	/**
	 * Returns all vertices in the underlying diagram graph, including site and
	 * constructed Voronoi vertices.
	 *
	 * @return every vertex currently stored in the diagram
	 */
	public List<Vertex> getVertices() {
		return new ArrayList<>(g.vertices);
	}

	/**
	 * Returns the face nearest to the given point.
	 *
	 * @param x x-coordinate of the query point
	 * @param y y-coordinate of the query point
	 * @return the nearest face
	 */
	public Face nearestFace(double x, double y) {
		return nearestFaces(x, y, 1).get(0);
	}

	/**
	 * Returns the N nearest faces to the given point.
	 *
	 * @param x x-coordinate of the query point
	 * @param y y-coordinate of the query point
	 * @param n maximum number of nearest faces to return; must be positive
	 * @return the nearest faces in ascending distance order
	 */
	public List<Face> nearestFaces(double x, double y, int n) {
		var heap = kdTree.findNearestNeighbors(new double[] { x, y }, n, DISTANCE_FUNCTION);
		var faces = new ArrayList<Face>(heap.size());

		while (heap.size() > 0) {
			faces.add(heap.getMax().face);
			heap.removeMax();
		}
		return faces;
	}

	/**
	 * Returns the number of user-inserted point sites in the diagram.
	 *
	 * @return the number of point sites, excluding the three bootstrap vertices
	 */
	public int numPointSites() {
		// the three initial vertices don't count
		return numPsites - 3;
	}

	/**
	 * Returns the number of inserted line sites.
	 *
	 * @return the number of line-segment sites
	 */
	public int numLineSites() {
		return numLsites;
	}

	/**
	 * Returns the number of inserted arc sites.
	 *
	 * @return the number of arc sites
	 */
	public int numArcSites() {
		return numAsites;
	}

	/**
	 * Returns the number of constructed Voronoi vertices.
	 * <p>
	 * Point-site vertices used as site handles are excluded from this count.
	 *
	 * @return the number of Voronoi vertices in the diagram
	 */
	public int numVertices() {
		return g.vertices.size() - numPointSites();
	}

	/**
	 * Returns the number of user-visible faces in the diagram.
	 *
	 * @return the number of non-null faces
	 */
	public int numFaces() {
		var count = 0;
		for (Face face : g.faces) {
			if (!face.isNullFace()) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Returns the total number of faces in the underlying graph.
	 *
	 * @return the number of all faces, including internal null faces
	 */
	public int numAllFaces() {
		return g.faces.size();
	}

	/**
	 * Returns the number of temporary split vertices currently present.
	 *
	 * @return the number of vertices with type {@link VertexType#SPLIT}
	 */
	public int numSplitVertices() {
		var count = 0;
		for (Vertex v : g.vertices) {
			if (v.type == VertexType.SPLIT) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Returns the underlying half-edge diagram used to store the Voronoi graph.
	 *
	 * @return the mutable backing half-edge diagram
	 */
	public HalfEdgeDiagram getDiagram() {
		return g;
	}

	/**
	 * Resets the global vertex numbering used by {@link Vertex}.
	 * <p>
	 * This is primarily useful in tests or tooling that expects deterministic
	 * vertex indices across independent diagram instances.
	 */
	public static void resetVertexCount() {
		Vertex.resetCount();
	}

	/**
	 * Runs the internal topology and geometry consistency checks.
	 *
	 * @return {@code true} if the current diagram passes the checker
	 */
	public boolean check() {
		return vd_checker.isValid();
	}

	/**
	 * Applies an edge filter to the underlying graph.
	 * <p>
	 * Any edge for which the filter returns {@code false} is marked invalid until
	 * {@link #filterReset()} is called.
	 *
	 * @param flt the filter to apply to every edge
	 */
	public void filter(Filter flt) {
		flt.setGraph(g);
		for (Edge e : g.edges) {
			if (!flt.apply(e)) {
				e.valid = false;
			}
		}
	}

	/**
	 * Clears any active edge filtering and marks every edge as valid again.
	 */
	public void filterReset() {
		for (Edge e : g.edges) {
			e.valid = true;
		}
	}

	// comparison-predicate for VertexQueue
	/**
	 * In {@code augmentVertexSet()} we grow the delete-tree by processing
	 * vertices one by one from a priority queue.
	 * <p>
	 * This is the priority queue sort predicate. We handle vertices with a large
	 * {@code fabs(in_circle())} first, since we believe their predicate to be more
	 * reliable.
	 */
	class AbsComparison implements Comparator<Pair<Vertex, Double>> {
		@Override
		public int compare(Pair<Vertex, Double> lhs, Pair<Vertex, Double> rhs) {
			return -Double.compare(Math.abs(lhs.getSecond()), Math.abs(rhs.getSecond()));
		}
	}

	// data required for adding a new edge
	/**
	 * Stores information related to a new edge while it is being added.
	 */
	class EdgeData {
		/** edge prior to v1 */
		Edge v1_prv;
		/** NEW edge source */
		Vertex v1;
		/** edge following v1 */
		Edge v1_nxt;
		/** edge prior to v2 */
		Edge v2_prv;
		/** NEW edge target */
		Vertex v2;
		/** edge following v2 */
		Edge v2_nxt;
		/** face of v1 and v2 */
		Face f;

		@Override
		public String toString() {
			return String.format("EdgeData(\n  v1_prv: %s\n  v1: %s\n  v1_nxt: %s\n  v2_prv: %s\n  v2: %s\n  v2_nxt: %s\n)", v1_prv, v1, v1_nxt, v2_prv, v2,
					v2_nxt);
		}
	};

	/**
	 * Initializes the diagram with three generators. Add one vertex at origin and
	 * three vertices at 'infinity' and their associated edges
	 */
	private void initialize() {
		var far_multiplier = 6D;
		// initial generators/sites:
		var gen1 = new Point(0, 3.0 * farRadius);
		var gen2 = new Point(-3.0 * Math.sqrt(3.0) * farRadius / 2.0, -3.0 * farRadius / 2.0);
		var gen3 = new Point(+3.0 * Math.sqrt(3.0) * farRadius / 2.0, -3.0 * farRadius / 2.0);
		// initial vd-vertices
		var vd1 = new Point(0, -3.0 * farRadius * far_multiplier);
		var vd2 = new Point(+3.0 * Math.sqrt(3.0) * farRadius * far_multiplier / 2.0, +3.0 * farRadius * far_multiplier / 2.0);
		var vd3 = new Point(-3.0 * Math.sqrt(3.0) * farRadius * far_multiplier / 2.0, +3.0 * farRadius * far_multiplier / 2.0);
		// add init vertices
		var v00 = g.addVertex(new Vertex(new Point(0, 0), VertexStatus.UNDECIDED, VertexType.NORMAL, gen1));
		var v01 = g.addVertex(new Vertex(vd1, VertexStatus.OUT, VertexType.OUTER, gen3));
		var v02 = g.addVertex(new Vertex(vd2, VertexStatus.OUT, VertexType.OUTER, gen1));
		var v03 = g.addVertex(new Vertex(vd3, VertexStatus.OUT, VertexType.OUTER, gen2));
		// add initial sites to graph
		var vert1 = g.addVertex(new Vertex(gen1, VertexStatus.OUT, VertexType.POINTSITE));
		var vert2 = g.addVertex(new Vertex(gen2, VertexStatus.OUT, VertexType.POINTSITE));
		var vert3 = g.addVertex(new Vertex(gen3, VertexStatus.OUT, VertexType.POINTSITE));

		// apex-points on the three edges:
		var a1 = g.addVertex(new Vertex(gen2.add(gen3).mult(0.5), VertexStatus.UNDECIDED, VertexType.APEX, gen2));
		var a2 = g.addVertex(new Vertex(gen1.add(gen3).mult(0.5), VertexStatus.UNDECIDED, VertexType.APEX, gen3));
		var a3 = g.addVertex(new Vertex(gen1.add(gen2).mult(0.5), VertexStatus.UNDECIDED, VertexType.APEX, gen1));

		// add face 1: v0-v1-v2 which encloses gen3
		var e1_1 = g.add_edge(v00, a1);
		var e1_2 = g.add_edge(a1, v01);
		var e2 = g.add_edge(v01, v02);
		var e3_1 = g.add_edge(v02, a2);
		var e3_2 = g.add_edge(a2, v00);
		var f1 = g.addFace();
		f1.setSite(new PointSite(gen3, f1, vert3));
		f1.setStatus(FaceStatus.NONINCIDENT);
		kdTree.addPoint(new double[] { gen3.x, gen3.y }, (PointSite) f1.getSite());
		g.setNextCycle(Arrays.asList(e1_1, e1_2, e2, e3_1, e3_2), f1, 1);

		// add face 2: v0-v02-v03 which encloses gen1
		var e4_1 = g.add_edge(v00, a2);
		var e4_2 = g.add_edge(a2, v02);
		var e5 = g.add_edge(v02, v03);
		var e6_1 = g.add_edge(v03, a3);
		var e6_2 = g.add_edge(a3, v00);
		var f2 = g.addFace();
		f2.setSite(new PointSite(gen1, f2, vert1));
		f2.setStatus(FaceStatus.NONINCIDENT);
		kdTree.addPoint(new double[] { gen1.x, gen1.y }, (PointSite) f2.getSite());
		g.setNextCycle(Arrays.asList(e4_1, e4_2, e5, e6_1, e6_2), f2, 1);

		// add face 3: v0-v3-v1 which encloses gen2
		var e7_1 = g.add_edge(v00, a3);
		var e7_2 = g.add_edge(a3, v03);
		var e8 = g.add_edge(v03, v01);
		var e9_1 = g.add_edge(v01, a1);
		var e9_2 = g.add_edge(a1, v00);
		var f3 = g.addFace();
		f3.setSite(new PointSite(gen2, f3, vert2)); // this constructor needs f3...
		f3.setStatus(FaceStatus.NONINCIDENT);
		kdTree.addPoint(new double[] { gen2.x, gen2.y }, (PointSite) f3.getSite());
		g.setNextCycle(Arrays.asList(e7_1, e7_2, e8, e9_1, e9_2), f3, 1);

		// set type.
		e1_1.type = EdgeType.LINE;
		e1_1.setParameters(f1.getSite(), f3.getSite(), false);
		e1_2.type = EdgeType.LINE;
		e1_2.setParameters(f1.getSite(), f3.getSite(), true);
		e2.type = EdgeType.OUTEDGE;
		e3_1.type = EdgeType.LINE;
		e3_1.setParameters(f2.getSite(), f1.getSite(), true);
		e3_2.type = EdgeType.LINE;
		e3_2.setParameters(f2.getSite(), f1.getSite(), false);
		e4_1.type = EdgeType.LINE;
		e4_1.setParameters(f2.getSite(), f1.getSite(), false);
		e4_2.type = EdgeType.LINE;
		e4_2.setParameters(f2.getSite(), f1.getSite(), true);
		e5.type = EdgeType.OUTEDGE;
		e6_1.type = EdgeType.LINE;
		e6_1.setParameters(f2.getSite(), f3.getSite(), false);
		e6_2.type = EdgeType.LINE;
		e6_2.setParameters(f2.getSite(), f3.getSite(), true);
		e7_1.type = EdgeType.LINE;
		e7_1.setParameters(f2.getSite(), f3.getSite(), true);
		e7_2.type = EdgeType.LINE;
		e7_2.setParameters(f2.getSite(), f3.getSite(), false);
		e8.type = EdgeType.OUTEDGE;
		e9_1.type = EdgeType.LINE;
		e9_1.setParameters(f1.getSite(), f3.getSite(), true);
		e9_2.type = EdgeType.LINE;
		e9_2.setParameters(f1.getSite(), f3.getSite(), false);

		// twin edges
		g.twinEdges(e1_1, e9_2);
		g.twinEdges(e1_2, e9_1);
		e2.twin = null; // the outermost edges have invalid twins
		e5.twin = null;
		e8.twin = null;
		g.twinEdges(e3_1, e4_2);
		g.twinEdges(e3_2, e4_1);
		g.twinEdges(e6_1, e7_2);
		g.twinEdges(e6_2, e7_1);

		assert (vd_checker.isValid()) : "diagram validity check failed";

	}

	/**
	 * find amount of clearance-disk violation on all vertices of face f
	 * 
	 * @return vertex with the largest clearance-disk violation
	 */
	private Vertex findSeedVertex(Face f, Site site) {
		var minPred = 0.0;
		Vertex minimalVertex = null;
		var first = true;
		var current = f.getEdge();
		var start = current;
		do {
			var q = current.target;
			if ((q.status != VertexStatus.OUT) && (q.type == VertexType.NORMAL)) {
				var h = q.inCircle(site.apexPoint(q.position));
				if (first || ((h < minPred) && (site.inRegion(q.position)))) {
					minPred = h;
					minimalVertex = q;
					first = false;
				}
			}
			current = current.next;
		} while (!current.equals(start));
		assert (minPred < 0) : "minPred < 0";
		return minimalVertex;
	}

	/**
	 * Given the set {@code inVertices} of IN vertices, find and return the adjacent
	 * IN-OUT edges.
	 * <p>
	 * Later NEW vertices are inserted into each of the found IN-OUT edges.
	 */
	private List<Edge> findInOutEdges() {
		assert (!inVertices.isEmpty()) : "inVertices must not be empty when searching for IN-OUT edges";
		List<Edge> output = new ArrayList<>(); // new vertices generated on these edges
		for (Vertex v : inVertices) {
			assert (v.status == VertexStatus.IN) : "v.status == VertexStatus.IN";
			for (Edge e : v.outEdges) {
				if (e.target.status == VertexStatus.OUT) {
					output.add(e); // this is an IN-OUT edge
				}
			}
		}
		assert (!output.isEmpty()) : "no IN-OUT edges found despite non-empty inVertices";
		return output;
	}

	// find EdgeData for a new edge
	/**
	 * On a face which has IN and OUT-vertices, find the sequence
	 * OUT-OUT-OUT-..-OUT-NEW(v1)-IN-...-IN-NEW(v2)-OUT-OUT and return
	 * {@code v1}/{@code v2} together with their previous and next edges.
	 *
	 * @param f face on which we search for vertices
	 * @param startverts contains NEW-vertices already found, which are not valid
	 *                   for this call to findEdgeData
	 * @param segment contains ENDPOINT vertices when we are inserting a
	 *                line-segment; these vertices are needed to ensure finding
	 *                correct points around sites/null-edges
	 */
	private EdgeData findEdgeData(Face f, List<Vertex> startverts, Pair<Vertex, Vertex> segment) {
		var ed = new EdgeData();
		ed.f = f;
		var current_edge = f.getEdge(); // start on some edge of the face
		var start_edge = current_edge;
		var found = false;
		do { // find OUT-NEW-IN vertices in this loop
			var next_edge = current_edge.next;

			var previous_vertex = current_edge.source;
			var current_vertex = current_edge.target;
			var next_vertex = next_edge.target;

			if (isV1Candidate(current_vertex, previous_vertex, next_vertex, segment)
					&& !startverts.contains(current_vertex)) {
				ed.v1 = current_vertex;
				ed.v1_prv = current_edge;
				ed.v1_nxt = current_edge.next;
				found = true;
			}
			current_edge = current_edge.next;
		} while (!current_edge.equals(start_edge) && !found);
		assert (found) : "v1 (OUT-NEW-IN) not found";

		// now search for v2
		start_edge = current_edge; // note that this search starts where we ended in the loop above!
		found = false;
		do { // find IN-NEW-OUT vertices in this loop
			var current_vertex = current_edge.target;
			if (current_vertex.status == VertexStatus.NEW && current_vertex.type != VertexType.SEPPOINT) {
				if (!current_vertex.equals(ed.v1)) { // -IN-NEW(v2)
					ed.v2 = current_vertex;
					ed.v2_prv = current_edge;
					ed.v2_nxt = current_edge.next;
					found = true;
				}
			}
			current_edge = current_edge.next;
		} while (!current_edge.equals(start_edge) && !found);
		assert (found) : "v2 (IN-NEW-OUT) not found";
		return ed;
	}

	/**
	 * Tests whether a vertex is a valid candidate for the v1 (OUT-NEW-IN) position
	 * in the edge-data search.  A vertex qualifies if it is NEW (not a SEPPOINT)
	 * and either the previous vertex is OUT/UNDECIDED and not a segment endpoint,
	 * or the next vertex is a segment endpoint.
	 */
	private boolean isV1Candidate(Vertex current, Vertex previous, Vertex next, Pair<Vertex, Vertex> segment) {
		if (current.status != VertexStatus.NEW || current.type == VertexType.SEPPOINT) {
			return false;
		}
		var previousNotEndpoint = !previous.equals(segment.getFirst()) && !previous.equals(segment.getSecond());
		var previousIsOutOrUndecided = previous.status == VertexStatus.OUT || previous.status == VertexStatus.UNDECIDED;
		var nextIsEndpoint = next.equals(segment.getFirst()) || next.equals(segment.getSecond());
		return (previousIsOutOrUndecided && previousNotEndpoint) || nextIsEndpoint;
	}

	// find and return edges on which we potentially need SPLIT vertices
	/**
	 * Walk around the face {@code f} and return edges whose endpoints are on
	 * separate sides of the {@code pt1-pt2} line.
	 * <p>
	 * TODO: not all edges found like this need SPLIT vertices, but it does not
	 * hurt to insert SPLIT-vertices in this case.
	 */
	private List<Edge> findSplitEdges(Face f, Point pt1, Point pt2) {
		assert (checkFace(f)) : "face check failed";
		List<Edge> out = new ArrayList<>();
		var current_edge = f.getEdge();
		var start_edge = current_edge;
		// int count=0;
		do { // FIND ALL! not just one.
			var trg = current_edge.target;
			var src = current_edge.source;
			var src_is_right = src.position.is_right(pt1, pt2);
			var trg_is_right = trg.position.is_right(pt1, pt2);
			if (src.type == VertexType.NORMAL || src.type == VertexType.APEX || src.type == VertexType.SPLIT) { // ?
																												// check
																												// edge-type
																												// instead?
				if (src_is_right != trg_is_right) {
					out.add(current_edge);
					assert (checkEdge(current_edge)) : "edge check failed during split-edge search";
				}
			}
			current_edge = current_edge.next;
		} while (!current_edge.equals(start_edge));

		return out;
	}

	// find a SPLIT vertex on the Face f
	// return true, and set v, if found.
	private Vertex findSplitVertex(Face f) {
		for (Vertex q : g.faceVertices(f)) {
			if (q.type == VertexType.SPLIT) {
				return q;
			}
		}
		return null;
	}

	// find an adjacent vertex to v, along an edge that is not a NULLEDGE
	// return true if found, otherwise false.
	private Vertex nullVertexTarget(Vertex v) {
		for (Edge e : v.outEdges) {
			if (e.type != EdgeType.NULLEDGE) {
				return e.target;
			}
		}
		return null;
	}

	/**
	 * Grows the delete-tree of IN vertices by "weighted breadth-first search"
	 * <p>
	 * We start at the seed and add vertices with {@code detH < 0}, provided that:
	 * <ul>
	 * <li>(C4) {@code v} should not be adjacent to two or more IN vertices, since
	 * that would result in a loop/cycle</li>
	 * <li>(C5) for an incident face containing {@code v}, {@code v} is adjacent to
	 * an IN vertex on this face</li>
	 * </ul>
	 * <p>
	 * C4 and C5 refer to the Sugihara&amp;Iri 1992 "one million" paper. We process
	 * UNDECIDED vertices adjacent to known IN-vertices in a "weighted
	 * breadth-first-search" manner where vertices with a large {@code fabs(detH)}
	 * are
	 * processed first, since we assume the in-circle predicate to be more reliable
	 * the larger {@code fabs(in_circle())} is.
	 * 
	 * @param site site currently being inserted
	 */
	private void augmentVertexSet(Site site) {
		while (!vertexQueue.isEmpty()) {
			var nextVertexDet = vertexQueue.poll();
			var v = nextVertexDet.getFirst();
			double h = nextVertexDet.getSecond();
		assert (v.status == VertexStatus.UNDECIDED) : "v.status == VertexStatus.UNDECIDED";
			if (h < 0.0) { // try to mark IN if h<0 and passes (C4) and (C5) tests and in_region().
							// otherwise mark OUT
				if (predicateC4(v) || !predicateC5(v) || !site.inRegion(v.position)) {
					v.status = VertexStatus.OUT; // C4 or C5 violated, so mark OUT
				} else {
					markVertex(v, site); // h<0 and no violations, so mark IN. push adjacent UNDECIDED vertices onto Q.
				}
			} else {
				v.status = VertexStatus.OUT; // detH was positive (or zero), so mark OUT
			}
			modifiedVertices.add(v);
		}

		assert (vertexQueue.isEmpty()) : "vertex queue should be empty after augmentation";
		assert (verticesAllIN(inVertices)) : "all inVertices should have IN status after augmentation";
		// sanity-check?: for all incident faces the IN/OUT-vertices should be connected
	}

	// adjacent in-count predicate for buildingdelete-tree
	/**
	 * Number of IN vertices adjacent to the given vertex {@code v}.
	 * <p>
	 * Predicate C4, i.e. "adjacent in-count", from the Sugihara&amp;Iri 1992
	 * "one million" paper.
	 */
	private boolean predicateC4(Vertex v) {
		var in_count = 0;
		for (Edge e : v.outEdges) {
			var w = e.target;
			if (w.status == VertexStatus.IN) {
				in_count++;
				if (in_count >= 2) {
					return true;
				}
			}
		}
		return false;
	}

	// connectedness-predicate for delete-tree building
	/**
	 * Do any of the three faces that are adjacent to the given IN-vertex
	 * {@code v} have an IN-vertex?
	 * <p>
	 * Predicate C5, i.e. "connectedness", from the Sugihara&amp;Iri 1992
	 * "one million" paper.
	 */
	private boolean predicateC5(Vertex v) {
		if (v.type == VertexType.APEX || v.type == VertexType.SPLIT) {
			return true;
		}

		Face f0 = null, f1 = null, f2 = null;
		int count = 0;

		for (Edge e : v.outEdges) {
			Face f = e.face;
			if (f.getStatus() != FaceStatus.INCIDENT) {
				continue;
			}
			if (f != f0 && f != f1 && f != f2) {
				if (count == 0) {
					f0 = f;
				} else if (count == 1) {
					f1 = f;
				} else {
					f2 = f;
				}
				count++;
			}
		}

		assert (count > 0) : "vertex must have at least one incident face";

		if (count > 0 && !faceHasAllowedNeighbor(v, f0)) {
			return false;
		}
		if (count > 1 && !faceHasAllowedNeighbor(v, f1)) {
			return false;
		}
		if (count > 2 && !faceHasAllowedNeighbor(v, f2)) {
			return false;
		}

		return true;
	}

	private boolean faceHasAllowedNeighbor(Vertex v, Face f) {
	    assert (v != null) : "vertex must not be null";
	    assert (f != null) : "face must not be null";

	    Edge out = null;
	    for (Edge e : v.outEdges) {
	        if (e.face == f) {
	            out = e;
	            break;
	        }
	    }

	    assert (out != null) : "v must have an outgoing edge on face f";
	    assert (out.prev != null) : "prev pointers must be maintained";
	    assert (out.prev.target == v) : "out.prev should enter v on face f";

	    Vertex next = out.target;
	    Vertex prev = out.prev.source;

	    if (isAllowedC5Neighbor(prev, v) || isAllowedC5Neighbor(next, v)) {
	        return true;
	    }

	    // Conservative fallback to preserve legacy behavior for special vertices
	    Edge current = out.next;
	    while (current != out) {
	        Vertex w = current.target;
	        if (w != v && (w.type == VertexType.ENDPOINT
	                || w.type == VertexType.APEX
	                || w.type == VertexType.SPLIT)) {
	            return true;
	        }
	        current = current.next;
	    }

	    return false;
	}

	private boolean isAllowedC5Neighbor(Vertex w, Vertex v) {
		if (w == null || w == v) {
			return false;
		}

		if (w.status == VertexStatus.IN) {
			return true;
		}

		if (w.type == VertexType.ENDPOINT || w.type == VertexType.APEX || w.type == VertexType.SPLIT || w.type == VertexType.SEPPOINT) {
			return true;
		}

		return false;
	}

	// mark adjacent faces INCIDENT
	// call this when inserting line-sites
	// since we call addSplitVertex we can't use iterators, because they get
	// invalidated
	// so use the slower adjacent_faces() instead.
	private void markAdjacentFaces(Vertex v, Site site) {
		assert (v.status == VertexStatus.IN) : "vertex must have IN status";
		var new_adjacent_faces = g.adjacentFaces(v);

		assert ((v.type == VertexType.APEX && new_adjacent_faces.size() == 2)
				|| (v.type == VertexType.SPLIT && new_adjacent_faces.size() == 2)
				|| new_adjacent_faces.size() == 3) : "unexpected adjacent face count for vertex type";

		for (Face adj_face : new_adjacent_faces) {
			if (adj_face.getStatus() != FaceStatus.INCIDENT) {
				if (site.isLine()) {
					addSplitVertex(adj_face, site);
				}

				adj_face.setStatus(FaceStatus.INCIDENT);
				incidentFaces.add(adj_face);
			}
		}
	}

	// mark adjacent faces INCIDENT
	// IN-Vertex v has three adjacent faces, mark nonincident faces incident
	// and push them to the incident_faces queue
	// NOTE: call this only when inserting point-sites
	private void markAdjacentFacesP(Vertex v) {
		assert (v.status == VertexStatus.IN) : "vertex must have IN status";
		for (Edge e : v.outEdges) {
			var adj_face = e.face;
			if (adj_face.getStatus() != FaceStatus.INCIDENT) {
				adj_face.setStatus(FaceStatus.INCIDENT);
				incidentFaces.add(adj_face);
			}
		}
	}

	/**
	 * mark vertex IN and mark adjacent faces INCIDENT; push adjacent UNDECIDED
	 * vertices onto queue
	 */
	private void markVertex(Vertex v, Site site) {
		v.status = VertexStatus.IN;
		inVertices.add(v);
		modifiedVertices.add(v);

		if (site.isPoint()) {
			markAdjacentFacesP(v);
		} else {
			markAdjacentFaces(v, site);
		}

		// push the v-adjacent vertices onto the queue
		for (Edge e : v.outEdges) {
			var w = e.target;
			if ((w.status == VertexStatus.UNDECIDED) && (!w.inQueue)) {
				// when pushing onto queue we also evaluate in_circle predicate so that we
				// process vertices in the correct order
				vertexQueue.add(new Pair<Vertex, Double>(w, w.inCircle(site.apexPoint(w.position))));
				w.inQueue = true;
			}
		}
	}

	/**
	 * Add NEW vertices on IN-OUT edges.
	 * <p>
	 * This generates new Voronoi vertices on all IN-OUT edges.
	 * 
	 * @implNote used only by {@link #insertPointSite(Point)}
	 * @param new_site site currently being inserted
	 */
	private void addVertices(Site new_site) {
		assert (!inVertices.isEmpty()) : "inVertices must not be empty when adding new vertices";
		var q_edges = findInOutEdges(); // new vertices generated on these IN-OUT edges
		for (Edge e : q_edges) {
			var sl = vpos.position(e, new_site); // vertex_positioner.cpp
			var q = g.addVertex(new Vertex(sl.p, VertexStatus.NEW, VertexType.NORMAL, new_site.apexPoint(sl.p), sl.k3));
			modifiedVertices.add(q);
			g.addVertexInEdge(q, e);
			q.maxError = vpos.distError(e, sl, new_site);
		}
	}

	// add a new face corresponding to the new Site
	/**
	 * Call add_new_edge() on all the incident_faces that should be split.
	 */
	private Face addFace(Site s) {
		var newface = g.addFace();
		newface.setSite(s);
		s.face = newface;
		newface.setStatus(FaceStatus.NONINCIDENT);
		if (s.isPoint()) {
			kdTree.addPoint(new double[] { s.position().x, s.position().y }, (PointSite) s);
		}
		return newface;
	}

	private void addEdges(Face newface, Face f) {
		addEdges(newface, f, null, new Pair<Vertex, Vertex>(null, null));
	}

	// add all NEW-NEW edges
	/**
	 * By adding a NEW-NEW edge, split the face {@code f} into one part which is
	 * {@code newface}, and the other part is the old {@code f}.
	 * <p>
	 * For linesegment or arc sites we pass in both the {@code k=+1} face
	 * {@code newface} and the {@code k=-1} face {@code newface2}. The segment
	 * endpoints are passed to {@link #findEdgeData}.
	 */
	private void addEdges(Face newface, Face f, Face newface2, Pair<Vertex, Vertex> segment) {
		var new_count = numNewVertices(f);
		assert (new_count > 0) : "new_count > 0";
		assert ((new_count % 2) == 0) : "(new_count % 2) == 0";
		var new_pairs = new_count / 2; // we add one NEW-NEW edge for each pair found
		List<Vertex> startverts = new ArrayList<>(); // this holds ed.v1 vertices for edges already added
		for (var m = 0; m < new_pairs; m++) {
			var ed = findEdgeData(f, startverts, segment);
			addEdge(ed, newface, newface2);
			startverts.add(ed.v1);
		}
	}

	// add a new edge to the diagram
	// newface = the k=+1 positive offset face
	// newface2 = the k=-1 negative offset face
	private void addEdge(EdgeData ed, Face newface, Face newface2) {
		var new_previous = ed.v1_prv;
		var new_source = ed.v1; // -OUT-NEW(v1)-IN-...
		var twin_next = ed.v1_nxt;

		var twin_previous = ed.v2_prv;
		var new_target = ed.v2; // -IN-NEW(v2)-OUT-
		var new_next = ed.v2_nxt;

		var f = ed.f;
		var f_site = f.getSite();
		Site new_site;
		Face new_face;
		if (new_source.k3 == 1) { // find out if newface or newface2 should be used
			new_site = newface.getSite();
			new_face = newface;
		} else {
			new_site = newface2.getSite();
			new_face = newface2;
		}

		// both trg and src should be on same side of new site
		assert (new_target.k3 == new_source.k3) : "new_target.k3 == new_source.k3";

		// Determine which side of the bisector the source and target vertices are on.
		// If they are on different sides, an apex-split vertex is required.
		var signs = determineApexSigns(new_source, new_target, f_site, new_site);
		var src_sign = signs[0];
		var trg_sign = signs[1];

		if (src_sign == trg_sign) {
			addSingleEdge(new_source, new_target, new_previous, new_next, twin_previous, twin_next,
					f, f_site, new_site, new_face, src_sign);
		} else {
			addApexSplitEdges(new_source, new_target, new_previous, new_next, twin_previous, twin_next,
					f, f_site, new_site, new_face, src_sign, trg_sign);
		}
	}

	/**
	 * Determines whether the source and target vertices of a new edge are on the
	 * same side ("sign") of the bisector between the face site and the new site.
	 * <p>
	 * When signs differ, an apex-split is required to insert two edges instead of one.
	 *
	 * @return a boolean array {@code [src_sign, trg_sign]}
	 */
	private boolean[] determineApexSigns(Vertex new_source, Vertex new_target, Site f_site, Site new_site) {
		var src_sign = true;
		var trg_sign = true;
		if (f_site.isPoint() && new_site.isLine()) { // PL
			var pt2 = f_site.position();
			var pt1 = new_site.apexPoint(pt2); // projection of pt2 onto LineSite
			src_sign = new_source.position.is_right(pt1, pt2);
			trg_sign = new_target.position.is_right(pt1, pt2);
		} else if (f_site.isPoint() && new_site.isArc()) { // PA
			var pt2 = f_site.position();
			var cen = new Point(new_site.x(), new_site.y());
			var cen_pt2 = pt2.sub(cen);
			var pt1 = cen.add(cen_pt2.mult(new_site.r() / cen_pt2.norm()));
			src_sign = new_source.position.is_right(pt1, pt2);
			trg_sign = new_target.position.is_right(pt1, pt2);
		} else if (f_site.isPoint() && new_site.isPoint()) { // PP
			src_sign = new_source.position.is_right(f_site.position(), new_site.position());
			trg_sign = new_target.position.is_right(f_site.position(), new_site.position());
		} else if (f_site.isLine() && new_site.isLine()) { // LL
			// a line-line bisector, sign should not matter because there is no sqrt()
			assertLineLineInRegion(new_source, new_target, f_site, new_site);
		} else if (f_site.isLine() && new_site.isArc()) { // LA
			var pt2 = new Point(new_site.x(), new_site.y());
			var pt1 = f_site.apexPoint(pt2);
			src_sign = new_source.position.is_right(pt1, pt2);
			trg_sign = new_target.position.is_right(pt1, pt2);
			// if one vertex is on a null-face, we cannot trust the sign
			if (new_source.dist() == 0 || new_target.dist() == 0) {
				if (new_source.dist() > new_target.dist()) {
					src_sign = trg_sign;
				} else {
					trg_sign = src_sign;
				}
			}
		} else {
			throw new RuntimeException(String.format(
					"addEdge: unhandled site-type combination f_site=%s, new_site=%s", f_site, new_site));
		}
		return new boolean[] { src_sign, trg_sign };
	}

	/**
	 * Asserts that both edge endpoints lie on the correct side of the two line
	 * sites in a line-line bisector configuration.  Only runs the check when
	 * the vertices are sufficiently far from the line endpoints.
	 */
	private void assertLineLineInRegion(Vertex new_source, Vertex new_target, Site f_site, Site new_site) {
		if ((new_source.position != new_target.position)
				&& (new_source.position != f_site.start())
				&& (new_source.position != f_site.end())
				&& (new_target.position != f_site.start())
				&& (new_target.position != f_site.end())
				&& (new_source.position.sub(f_site.apexPoint(new_source.position)).norm() > Numeric.DISTANCE_EPSILON)
				&& (new_target.position.sub(f_site.apexPoint(new_target.position)).norm() > Numeric.DISTANCE_EPSILON)) {
			assert (!new_source.position.is_right(f_site.start(), f_site.end())) : "!new_source.position.is_right(f_site)";
			assert (!new_target.position.is_right(f_site.start(), f_site.end())) : "!new_target.position.is_right(f_site)";
			assert (!new_source.position.is_right(new_site.start(), new_site.end())) : "!new_source.position.is_right(new_site)";
			assert (!new_target.position.is_right(new_site.start(), new_site.end())) : "!new_target.position.is_right(new_site)";
		}
	}

	/**
	 * Adds a single NEW-NEW edge (no apex-split required because both vertices
	 * are on the same side of the bisector).
	 */
	private void addSingleEdge(Vertex new_source, Vertex new_target,
			Edge new_previous, Edge new_next, Edge twin_previous, Edge twin_next,
			Face f, Site f_site, Site new_site, Face new_face, boolean src_sign) {
		var twinEdges = g.addTwinEdges(new_source, new_target);
		var e_new = twinEdges.getFirst();
		var e_twin = twinEdges.getSecond();
		g.setNext(e_new, new_next);
		assert (new_next.k == new_previous.k) : "new_next.k == new_previous.k";
		e_new.k = new_next.k; // the next edge is on the same face, so has the correct k-value
		e_new.face = f; // src-trg edge has f on its left
		g.setNext(new_previous, e_new);
		f.setEdge(e_new);
		e_new.setParameters(f_site, new_site, !src_sign);

		g.setNext(twin_previous, e_twin);
		g.setNext(e_twin, twin_next);
		e_twin.k = new_source.k3;
		e_twin.setParameters(f_site, new_site, !src_sign);
		e_twin.face = new_face;
		new_face.setEdge(e_twin);

		assert (checkEdge(e_new) && checkEdge(e_twin)) : "checkEdge(e_new) && checkEdge(e_twin)";
	}

	/**
	 * Adds two NEW-APEX-NEW edges when the source and target are on different
	 * sides of the bisector (apex-split is required).
	 */
	private void addApexSplitEdges(Vertex new_source, Vertex new_target,
			Edge new_previous, Edge new_next, Edge twin_previous, Edge twin_next,
			Face f, Site f_site, Site new_site, Face new_face, boolean src_sign, boolean trg_sign) {
		//     f                          f
		// new_prv -> NEW -- e1 ----> APEX --e2 ---> NEW -> new_nxt
		// twn_nxt <- NEW <- e1_tw -- APEX <-e2_tw-- NEW <- twn_prv
		//            new1/new2               new1/new2
		var apex = g.addVertex(new Vertex(new Point(0, 0), VertexStatus.NEW, VertexType.APEX));
		var twin_edges1 = g.addTwinEdges(new_source, apex);
		var e1 = twin_edges1.getFirst();
		var e1_tw = twin_edges1.getSecond();
		var twin_edges2 = g.addTwinEdges(apex, new_target);
		var e2 = twin_edges2.getFirst();
		var e2_tw = twin_edges2.getSecond();
		e1.setParameters(f_site, new_site, !src_sign);
		e2.setParameters(f_site, new_site, !trg_sign);

		assert (new_previous.face == f) : "new_previous.face == f";
		assert (new_next.face == new_previous.face) : "new_next.face == new_previous.face";
		assert (new_next.k == new_previous.k) : "new_next.k == new_previous.k";

		// new_previous -> e1 -> e2 -> new_next
		g.setNext(new_previous, e1);
		g.setNext(e1, e2);
		g.setNext(e2, new_next);
		e1.face = f;
		e2.face = f;
		e1.k = new_next.k;
		e2.k = new_next.k;
		f.setEdge(e1);
		// twin edges
		e1_tw.setParameters(new_site, f_site, src_sign);
		e2_tw.setParameters(new_site, f_site, trg_sign);

		assert (twin_previous.k == twin_next.k) : "twin_previous.k == twin_next.k";
		assert (twin_previous.face == twin_next.face) : "twin_previous.face == twin_next.face";
		// twin_prev -> e2_tw -> e1_tw -> twin_next on new_face

		g.setNext(twin_previous, e2_tw);
		g.setNext(e2_tw, e1_tw);
		g.setNext(e1_tw, twin_next);

		e1_tw.k = new_source.k3;
		e2_tw.k = new_source.k3;
		new_face.setEdge(e1_tw);
		e1_tw.face = new_face;
		e2_tw.face = new_face;

		assert (checkEdge(e1) && checkEdge(e1_tw)) : "checkEdge(e1) && checkEdge(e1_tw)";
		assert (checkEdge(e2) && checkEdge(e2_tw)) : "checkEdge(e2) && checkEdge(e2_tw)";

		// position the apex
		var min_t = e1.minimumT(f_site, new_site);
		apex.position = e1.point(min_t);
		apex.initDist(f_site.apexPoint(apex.position));
		modifiedVertices.add(apex);
	}

	class SeparatorTarget {
		public Edge vPrevious;
		public Vertex v_target;
		public Edge vNext;
		public boolean outNewIn;

		public SeparatorTarget(Edge vPrevious, Vertex v_target, Edge vNext, boolean outNewIn) {
			this.vPrevious = vPrevious;
			this.v_target = v_target;
			this.vNext = vNext;
			this.outNewIn = outNewIn;
		}

		@Override
		public String toString() {
			return String.format("SeparatorTarget(%s, %s, %s, %s)", vPrevious, v_target, vNext, outNewIn);
		}
	}

	/**
	 * Add SEPARATOR edge on the face f, which contains the endpoint
	 * 
	 * @param f         face of endpoint
	 * @param null_face null face of the endpoint
	 * @param target    target-data found by ??
	 * @param sep_endp
	 * @param s1        positive LineSite
	 * @param s2        negative LineSite
	 */
	private void addSeparator(Face f, Face nullFace, SeparatorTarget target, Vertex sep_endp, Site s1, Site s2) {
		if (sep_endp == null) {
			return; // do nothing!
		}

		assert ((sep_endp.k3 == 1) || (sep_endp.k3 == -1)) : "(sep_endp.k3 == 1) || (sep_endp.k3 == -1)";
		sep_endp.zeroDist();

		var next_prev = g.findNextPrev(nullFace, sep_endp);
		var endp_next_tw = next_prev.getFirst();
		var endp_prev_tw = next_prev.getSecond();
		var endp_prev = endp_next_tw.twin; // NOTE twin!
		var endp_next = endp_prev_tw.twin; // NOTE twin!
		assert (endp_next != null) : "endp_next != null";
		assert (endp_prev != null) : "endp_prev != null";

		// find NEW vertex on the old face f
		// this vertex has the correct alfa angle for this endp/separator
		var vPrevious = target.vPrevious;
		var v_target = target.v_target;
		var vNext = target.vNext;
		var outNewIn = target.outNewIn;
		assert ((v_target.k3 == 1) || (v_target.k3 == -1)) : "(v_target.k3 == 1) || (v_target.k3 == -1)";
		assert (sep_endp.k3 == v_target.k3) : "sep_endp.k3 == v_target.k3";
		// can't assert about in_region - numerical error is always present
		// assert( s1.in_region(v_target.position ) ) : " s1.in_region(v_target.position
		// ) ";
		// assert( s2.in_region(v_target.position ) ) : " s2.in_region(v_target.position
		// ) ";

		// add new separator edge, and its twin
		var twinEdges = g.addTwinEdges(sep_endp, v_target);
		var e2 = twinEdges.getFirst();
		var e2_tw = twinEdges.getSecond();
		e2.type = EdgeType.SEPARATOR;
		e2_tw.type = EdgeType.SEPARATOR;

		// there are two cases. depending on how v_target (NEW) is found:
		// OUT-NEW-IN, when out_new_in = true
		// IN-NEW-OUT, when out_new_in = false
		// here we set the faces, sites, and next-pointers depending on the case
		if (outNewIn) {
			e2.k = v_target.k3; // e2 is on the segment side
			e2_tw.k = +1; // e2_tw is on the point-site side

			e2_tw.face = f; // point-site face
			e2_tw.nullFace = f;
			e2_tw.hasNullFace = true;

			f.setEdge(e2_tw);
			endp_prev.k = e2.k; // endp_prev is on the line-site side

			if (e2.k == -1) { // choose either s1 or s2 as the site
				e2.face = s2.face;
				s2.face.setEdge(e2);
				endp_prev.face = s2.face;
			} else {
				e2.face = s1.face;
				s1.face.setEdge(e2);
				endp_prev.face = s1.face;
			}

			g.setNext(vPrevious, e2_tw);
			g.setNext(e2_tw, endp_next);

			endp_next.face = f; // the null-edge
			endp_next.k = 1;

			g.setNext(e2, vNext);
		} else {
			e2.k = +1; // e2 is on the point-site side
			e2_tw.k = v_target.k3; // e2_tw is on the segment side

			e2.face = f; // point-site face
			e2.nullFace = f;
			e2.hasNullFace = true;

			f.setEdge(e2);
			endp_next.k = e2_tw.k; // endp_next is on the linesite-side
			if (e2_tw.k == -1) {
				e2_tw.face = s2.face;
				s2.face.setEdge(e2_tw);
				endp_next.face = s2.face;
			} else {
				e2_tw.face = s1.face;
				s1.face.setEdge(e2_tw);
				endp_next.face = s1.face;
			}
			g.setNext(vPrevious, e2_tw);
			endp_prev.face = f;
			endp_prev.k = 1;

			g.setNext(endp_prev, e2);
			g.setNext(e2, vNext);
		}
		e2.setSepParameters(sep_endp.position, v_target.position);
		e2_tw.setSepParameters(sep_endp.position, v_target.position);

		assert (checkEdge(e2)) : "edge check failed for separator edge";
		assert (checkEdge(e2_tw)) : "edge check failed for separator twin edge";

	}

	/**
	 * Adds one or many SPLIT vertices to the edges of the give face. These are
	 * projections/mirrors of the site of f with the new Site s acting as the
	 * mirror. SPLIT vertices are inserted to avoid deleting loops during //
	 * {@link #augmentVertexSet}.
	 * 
	 * @param f
	 * @param s
	 */
	private void addSplitVertex(Face f, Site s) {
		assert (!s.isPoint()) : "split vertices can only be added for line/arc sites, not point sites";

		var fs = f.getSite();

		// don't search for split-vertex on the start or end face
		if (fs.isPoint() && s.isLine()) {
			if (fs.position() == s.start() || fs.position() == s.end()) {
				return;
			}
		}

		if (fs.isPoint() && s.isLine() && s.inRegion(fs.position())) {
			// 1) find the correct edge
			var pt1 = fs.position();
			var pt2 = pt1.sub(new Point(s.a(), s.b()));

			assert ((pt1.sub(pt2)).norm() > 0) : "(pt1.sub(pt2)).norm() > 0";

			var split_edges = findSplitEdges(f, pt1, pt2);
			// the sought edge should have src on one side of pt1-pt2
			// and trg on the other side of pt1-pt2

			for (Edge split_edge : split_edges) {
				if ((split_edge.type == EdgeType.SEPARATOR) || (split_edge.type == EdgeType.LINESITE)) {
					return; // don't place split points on linesites or separators(?)
					// find a point = src + u*(trg-src)
					// with min_t < u < max_t
					// and minimum distance to the pt1-pt2 line
				}

				Point split_pt_pos;

				var split_src = split_edge.source;
				var split_trg = split_edge.target;
				var errFunctr = new SplitPointError(g, split_edge, pt1, pt2); // error functor
				var min_t = Math.min(split_src.dist(), split_trg.dist());
				var max_t = Math.max(split_src.dist(), split_trg.dist());
				// require that min_t and max_t bracket the root
				if (errFunctr.value(min_t) * errFunctr.value(max_t) >= 0) {
					return;
				}

				var solver = new BracketingNthOrderBrentSolver(Numeric.BRENT_SOLVER_ABSOLUTE_EPSILON, 5);
				var max_iter = 500;
				var result = solver.solve(max_iter, errFunctr, min_t, max_t, AllowedSolution.ANY_SIDE);

				split_pt_pos = split_edge.point(result);

				var v = g.addVertex(new Vertex(split_pt_pos, VertexStatus.UNDECIDED, VertexType.SPLIT, fs.position()));

				assert (checkEdge(split_edge)) : "edge check failed for split edge";
				// 3) insert new SPLIT vertex into the edge
				g.addVertexInEdge(v, split_edge);
			}
		}
	}

	/**
	 * Result of finding or creating a null-face at a segment endpoint.
	 */
	private class NullFaceResult {

		/** new ENDPOINT vertex for the segment start/end */
		public final Vertex segEndpoint;
		/** null-face at the endpoint (new or existing) */
		public final Face nullFace;
		/** positive SEPARATOR edge endpoint vertex (if a positive separator should be added) */
		public final Vertex posSepVertex;
		/** negative SEPARATOR edge endpoint vertex (if a negative separator should be added) */
		public final Vertex negSepVertex;
		/** face-to-null: if a PointSite face should disappear, we return it here */
		public final Face faceToNull;

		public NullFaceResult(Vertex segEndpoint, Face nullFace, Vertex posSepVertex, Vertex negSepVertex, Face faceToNull) {
			this.segEndpoint = segEndpoint;
			this.nullFace = nullFace;
			this.posSepVertex = posSepVertex;
			this.negSepVertex = negSepVertex;
			this.faceToNull = faceToNull;
		}
	}

	// either find an existing null-face, or create a new one.
	/**
	 * @param start the end of the segment for which we will find/create a null-face
	 * @param other the other end of the new segment
	 * @param left a point left of the new segment
	 * @param dir alfa-direction for positioning endpoint vertex on null-face
	 * @param new_site the new Site we are inserting
	 * Returns:
	 * <ul>
	 * <li>HEVertex ENDPOINT-vertex for the new vertex</li>
	 * <li>HEFace null-face at endpoint (new or existing)</li>
	 * <li>HEVertex positive SEPARATOR edge endpoint vertex (if a positive separator should be added)</li>
	 * <li>HEVertex negative SEPARATOR edge endpoint vertex (if a negative separator should be added)</li>
	 * <li>HEFace face-to-null. if a PointSite face should disappear, we return it here.</li>
	 * </ul>
	 */
	private NullFaceResult findNullFace(Vertex start, Vertex other, Point left, Point dir, Site new_site) {
		Vertex seg_start = null; // new end-point vertices
		Face start_null_face = null; // either existing or new null-face at start-vertex

		Vertex pos_sep_start = null; // optional separator endpoints at start
		Vertex neg_sep_start = null; // invalid vertices are default

		Face face_to_null = null; // invalid face is default

		// this works for LineSite
		var k3_sign = true;
		if (new_site.isLine()) {
			k3_sign = left.is_right(start.position, other.position);
		} else if (new_site.isArc()) {
			k3_sign = new Point(new_site.x(), new_site.y()).is_right(start.position, other.position);
		} else {
			throw new RuntimeException();
		}
		// this is used below and in findNullFace()
		// k3_sign is already calculated in insertLineSite() ??

		if (start.nullFace != null) { // there is an existing null face
			start_null_face = start.nullFace;

			// create a new segment ENDPOINT vertex with zero clearance-disk
			seg_start = g.addVertex(new Vertex(start.position, VertexStatus.OUT, VertexType.ENDPOINT, 0));
			// find the edge on the null-face where we insert seg_start
			Edge insert_edge = null;
			var current2 = start_null_face.getEdge();
			var start_edge2 = current2;
			seg_start.setAlfa(dir);
			var found = false;
			do {
				var face_incident = (current2.twin.face.getStatus() == FaceStatus.INCIDENT);
				if (face_incident) { // pick any incident face!
					insert_edge = current2;
					found = true;
				}
				current2 = current2.next;
			} while (!current2.equals(start_edge2) && !found); // FIXME end early with !found
			assert (insert_edge != null) : "insert_edge != null";
			assert (found) : "found";
			g.addVertexInEdge(seg_start, insert_edge); // insert endpoint in null-edge

			// "process" the adjacent null-edges
			var next_prev = g.findNextPrev(start_null_face, seg_start);
			var next_edge = next_prev.getFirst();
			var prev_edge = next_prev.getSecond();
			var neg_null_edge = processNullEdge(dir, next_edge, k3_sign, true);
			neg_sep_start = neg_null_edge.getFirst();
			face_to_null = neg_null_edge.getSecond();
			var pos_null_edge = processNullEdge(dir, prev_edge, k3_sign, false);
			pos_sep_start = pos_null_edge.getFirst();
			face_to_null = pos_null_edge.getSecond();
			return new NullFaceResult(seg_start, start_null_face, pos_sep_start, neg_sep_start, face_to_null);
		} else { // no existing null-face
			// create a new null face at start. the face has three vertices/edges:
			//
			// neg_sep -> seg_endp -> pos_sep
			//
			start_null_face = g.addFace(); // this face to the left of start->end edge
			start_null_face.setNullFace(true);

			seg_start = g.addVertex(new Vertex(start.position, VertexStatus.OUT, VertexType.ENDPOINT));
			seg_start.zeroDist();
			seg_start.setAlfa(dir);
			seg_start.k3 = 0;
			pos_sep_start = g.addVertex(new Vertex(start.position, VertexStatus.UNDECIDED, VertexType.SEPPOINT));
			neg_sep_start = g.addVertex(new Vertex(start.position, VertexStatus.UNDECIDED, VertexType.SEPPOINT));

			pos_sep_start.zeroDist();
			neg_sep_start.zeroDist();

			if (k3_sign) {
				pos_sep_start.k3 = +1;
				neg_sep_start.k3 = -1;
			} else {
				pos_sep_start.k3 = -1;
				neg_sep_start.k3 = +1;
			}
			pos_sep_start.setAlfa(dir.xyPerp().mult(+1));
			neg_sep_start.setAlfa(dir.xyPerp().mult(-1));
			// null-edges around the face
			var twin_edges1 = g.addTwinEdges(seg_start, pos_sep_start);
			var e1 = twin_edges1.getFirst();
			var e1_tw = twin_edges1.getSecond();
			var twin_edges2 = g.addTwinEdges(pos_sep_start, neg_sep_start);
			var e2 = twin_edges2.getFirst();
			var e2_tw = twin_edges2.getSecond();
			var twin_edges3 = g.addTwinEdges(neg_sep_start, seg_start);
			var e3 = twin_edges3.getFirst();
			var e3_tw = twin_edges3.getSecond();

			// e1 -> e2 -> e3 on start_null_face, k=1
			// e1t <- e2t <- e3t on g[start].face, k=1
			g.setNextCycle(Arrays.asList(e1, e2, e3), start_null_face, 1);
			var start_face = start.face;
			var start_face_edge = start_face.getEdge(); // crazy workaround, because set_next_cycles sets g[face].edge wrong
													// here!
			g.setNextCycle(Arrays.asList(e3_tw, e2_tw, e1_tw), start.face, 1);
			start_null_face.setEdge(e1);
			start_face.setEdge(start_face_edge);

			e1.type = EdgeType.NULLEDGE;
			e2.type = EdgeType.NULLEDGE;
			e3.type = EdgeType.NULLEDGE;
			e1_tw.type = EdgeType.NULLEDGE;
			e3_tw.type = EdgeType.NULLEDGE;
			e2_tw.type = EdgeType.NULLEDGE;

			start.nullFace = start_null_face;
			start_null_face.setSite(start_face.getSite());

			return new NullFaceResult(seg_start, start_null_face, pos_sep_start, neg_sep_start, face_to_null);
		}

	}

	/**
	 * Finds the target of a new SEPARATOR edge
	 * 
	 * @param f    the HEFace on which we search for the target vertex
	 * @param endp the end-point of the null-face with the SEPARATOR source
	 * @return
	 */
	private SeparatorTarget findSeparatorTarget(Face f, Vertex endp) {
		if (endp == null) {
			return new SeparatorTarget(null, null, null, false);
		}

		var current_edge = f.getEdge();
		var start_edge = current_edge;
		var found = false;
		Vertex v_target = null;
		Edge vPrevious = null, vNext = null;
		var outNewIn = true;
		do {
			var next_edge = current_edge.next;
			var previous_vertex = current_edge.source;
			var current_vertex = current_edge.target;
			var next_vertex = next_edge.target;
			var isOutNewIn = isOutNewInSequence(previous_vertex, current_vertex, next_vertex);
			var isInNewOut = isInNewOutSequence(previous_vertex, current_vertex, next_vertex);
			if ((isOutNewIn || isInNewOut) && (endp.k3 == current_vertex.k3) && !endp.equals(current_vertex)) {
				v_target = current_vertex;
				vPrevious = current_edge;
				vNext = next_edge;
				outNewIn = isOutNewIn;
				found = true;
			}
			current_edge = current_edge.next;
		} while (!current_edge.equals(start_edge) && !found);
		assert (found) : "separator target not found";

		return new SeparatorTarget(vPrevious, v_target, vNext, outNewIn);
	}

	/** Tests whether vertices follow the OUT/UNDECIDED - NEW - IN pattern. */
	private boolean isOutNewInSequence(Vertex prev, Vertex curr, Vertex next) {
		return (prev.status == VertexStatus.OUT || prev.status == VertexStatus.UNDECIDED)
				&& curr.status == VertexStatus.NEW
				&& next.status == VertexStatus.IN;
	}

	/** Tests whether vertices follow the IN - NEW - OUT/UNDECIDED pattern. */
	private boolean isInNewOutSequence(Vertex prev, Vertex curr, Vertex next) {
		return prev.status == VertexStatus.IN
				&& curr.status == VertexStatus.NEW
				&& (next.status == VertexStatus.OUT || next.status == VertexStatus.UNDECIDED);
	}

	// prepare null-face
	// next_edge lies on an existing null face
	// - Here we either insert a NEW NORMAL or SEPPOINT in the edge
	// - OR we push and convert an existing vertex.
	// the next_prev flag indicates if we are dealing with the next edge from the
	// new segment-endpoint next_prev=true
	// or if we are dealing with the previous edge (next_prev=false)
	private Pair<Vertex, Face> processNullEdge(Point dir, Edge next_edge, boolean k3, boolean next_prev) {
		assert (next_edge.type == EdgeType.NULLEDGE) : "next_edge.type == EdgeType.NULLEDGE";
		var trg = next_edge.target;
		var src = next_edge.source;

		var adj = next_prev ? trg : src; // this is the vertex adjacent to the end-point, on the null face
		assert ((next_prev ? src : trg).type == VertexType.ENDPOINT) : "(next_prev ? src : trg).type == VertexType.ENDPOINT";

		Vertex sep_point = null;
		double dir_mult = next_prev ? +1 : -1;
		var sep_dir = dir.xyPerp().mult(dir_mult);
		var sep_alfa = Numeric.diangle(sep_dir.x, sep_dir.y); // alfa of potential SEPPOINT

		double new_k3; // k3-value of a new or pushed vertex
		if (next_prev) {
			new_k3 = k3 ? +1 : -1;
		} else {
			new_k3 = k3 ? -1 : +1;
		}

		if (adj.type == VertexType.ENDPOINT) { // target is endpoint
			// insert a normal vertex, positioned at mid-alfa between src/trg.
			var new_v = g.addVertex(new Vertex(src.position, VertexStatus.NEW, VertexType.NORMAL, src.position));
			var mid = Numeric.diangleMid(src.alfa, trg.alfa);
			new_v.alfa = mid;
			modifiedVertices.add(new_v);
			g.addVertexInEdge(new_v, next_edge);
			new_v.k3 = new_k3;

			return new Pair<Vertex, Face>(null, null);
		} else {
			// Not an ENDPOINT vertex.
			assert (adj.type != VertexType.ENDPOINT) : "adjacent vertex on null-face should not be ENDPOINT here";
			var mid = 0D;
			var seppoint_pred = false;
			var parallel_pred = false;

			// set the two predicates
			if (next_prev) {
				var next_next = next_edge.next;
				var next_previous = g.previousEdge(next_edge);
				var next_trg = next_next.target; // source
				mid = Numeric.diangleMid(src.alfa, next_trg.alfa); // prev_src, trg
				seppoint_pred = (next_trg.type != VertexType.ENDPOINT);

				var next_out_trg = nullVertexTarget(next_edge.target);
				var prev_out_trg = nullVertexTarget(next_previous.source);
				if (next_out_trg != null && prev_out_trg != null) {
					parallel_pred = (((next_out_trg.status == VertexStatus.OUT) || (next_out_trg.status == VertexStatus.NEW)
							|| (next_out_trg.status == VertexStatus.UNDECIDED))
							&& ((prev_out_trg.status == VertexStatus.OUT) || (prev_out_trg.status == VertexStatus.NEW)
									|| (prev_out_trg.status == VertexStatus.UNDECIDED)));
				}
			} else { // !next_prev
				var prev_prev = g.previousEdge(next_edge);
				var next_next2 = next_edge.next;
				var prev_src = prev_prev.source;
				mid = Numeric.diangleMid(prev_src.alfa, trg.alfa);
				seppoint_pred = (prev_src.type != VertexType.ENDPOINT);

				var next_out_trg2 = nullVertexTarget(next_edge.source);
				var prev_out_trg2 = nullVertexTarget(next_next2.target);
				if (next_out_trg2 != null && prev_out_trg2 != null) {
					parallel_pred = (((next_out_trg2.status == VertexStatus.OUT) || (next_out_trg2.status == VertexStatus.NEW)
							|| (next_out_trg2.status == VertexStatus.UNDECIDED))
							&& ((prev_out_trg2.status == VertexStatus.OUT) || (prev_out_trg2.status == VertexStatus.NEW)
									|| (prev_out_trg2.status == VertexStatus.UNDECIDED)));
				}

			} // predicates now set.

			if (parallel_pred && adj.type == VertexType.SEPPOINT) {
				Edge sep_edge = null;

				for (Edge e : adj.outEdges) {
					assert (e.source == adj) : "outgoing edge source must match the adjacent vertex";
					if (e.type == EdgeType.SEPARATOR) {
						sep_edge = e;
					}
				}
				assert (sep_edge != null) : "separator edge must exist for SEPPOINT vertex";

				Edge sep_twin = sep_edge.twin;
				Face sep_face = sep_edge.face;
				var sep_site = sep_face.getSite();
				Face sep_twin_face = sep_twin.face;
				var sep_twin_site = sep_twin_face.getSite();

				Edge pointsite_edge = null;
				if (sep_site.isPoint()) {
					pointsite_edge = sep_edge;
				} else if (sep_twin_site.isPoint()) {
					pointsite_edge = sep_twin;
				}

				// set the separator target to NEW
				Vertex sep_target = sep_edge.target;
				sep_target.status = VertexStatus.NEW;
				sep_target.k3 = new_k3;
				modifiedVertices.add(sep_target);

				return new Pair<Vertex, Face>(null, pointsite_edge.face); // Returns the resolved face!
			}

			var adj_out = nullVertexTarget(adj);
			assert (adj_out != null) : "null-vertex target must exist for adjacent vertex";

			if (adj_out.status == VertexStatus.OUT || adj_out.status == VertexStatus.UNDECIDED) {
				sep_point = addSeparatorVertex(src, next_edge, sep_dir);
				sep_point.k3 = new_k3;
				return new Pair<Vertex, Face>(sep_point, null);
			} else {
				// target is not endpoint so we push and convert the vertex

				if (seppoint_pred) {
					// the pushed vertex becomes a SEPPOINT
					adj.alfa = sep_alfa;
					adj.type = VertexType.SEPPOINT;
					adj.status = VertexStatus.NEW;
					sep_point = adj;
				} else {
					// otherwise it becomes a normal NEW vertex
					adj.alfa = mid;
					adj.type = VertexType.NORMAL;
					adj.status = VertexStatus.NEW;
				}
				adj.k3 = new_k3;
				modifiedVertices.add(adj);
				return new Pair<Vertex, Face>(sep_point, null);
			}
		}
	}

	// add a SEPPOINT vertex into the given null-edge
	/**
	 * @param endp the endpoint corresponding to the null-edge/null-face
	 * @param edge the null-edge into which we insert the new vertex
	 * @param sep_dir direction for setting alfa of the new vertex
	 */
	private Vertex addSeparatorVertex(Vertex endp, Edge edge, Point sep_dir) {
		var sep = g.addVertex(new Vertex(endp.position, VertexStatus.OUT, VertexType.SEPPOINT));
		sep.setAlfa(sep_dir);
		g.addVertexInEdge(sep, edge);
		modifiedVertices.add(sep);
		return sep;
	}

	/**
	 * One-argument version of {@link #repairFace(Face, Pair, Pair, Pair)} used by
	 * {@link #insertPointSite(Point)}.
	 */
	private void repairFace(Face f) {
		repairFace(f, new Pair<Vertex, Vertex>(null, null), new Pair<Face, Face>(null, null), new Pair<Face, Face>(null, null));
	}

	/**
	 * Repairs next-pointers of HEFace
	 * 
	 * @param f
	 * @param segment
	 * @param nulled_faces
	 * @param null_face
	 */
	private void repairFace(Face f, Pair<Vertex, Vertex> segment, Pair<Face, Face> nulled_faces, Pair<Face, Face> nullFace) {
		// Walk around the face and repair the next-pointers.
		// Called on the newly created face after all NEW-NEW edges have been added.
		var current_edge = f.getEdge();
		var start_edge = current_edge;
		var c = 0;
		do {
			assert (checkEdge(current_edge)) : "checkEdge(current_edge)";
			var current_target = current_edge.target;
			var current_source = current_edge.source;
			for (Edge e : current_target.outEdges) {
				var out_target = e.target;
				if (!isRepairCandidate(out_target, current_source)) {
					continue;
				}

				// Override face-assignment for edges that need re-parenting
				if (shouldOverrideFaceAssignment(e, current_edge, current_source, current_target, out_target,
						segment, nulled_faces, nullFace)) {
					e.face = f;
					e.k = current_edge.k;
				}

				if (e.face == f) {
					g.setNext(current_edge, e);
					assert (current_edge.k == e.k) : "current_edge.k == e.k";
					assert (current_edge.face == current_edge.next.face) : "current_edge.face == current_edge.next.face";
				}
			}
			current_edge = current_edge.next; // jump to the next edge
			c++;
			if (c > 30000) {
				throw new AssertionError("face repair exceeded 30000 iterations — probable cycle in edge graph");
			}
		} while (!current_edge.equals(start_edge));
	}

	/**
	 * Tests whether an edge target is a valid candidate during face repair:
	 * it must not loop back to where we came from, and must be a NEW,
	 * ENDPOINT, or SEPPOINT vertex.
	 */
	private boolean isRepairCandidate(Vertex out_target, Vertex current_source) {
		if (out_target.equals(current_source)) {
			return false;
		}
		return out_target.status == VertexStatus.NEW
				|| out_target.type == VertexType.ENDPOINT
				|| out_target.type == VertexType.SEPPOINT;
	}

	/**
	 * Determines whether the face-assignment of edge {@code e} should be
	 * overridden during face repair.  This happens in two cases:
	 * <ol>
	 * <li>The edge is a null-edge entering/leaving separator or endpoint vertices
	 *     that is not on a null-face.</li>
	 * <li>The edge belongs to a face that has been "nulled" (i.e. a point-site
	 *     face that disappeared during line-segment insertion).</li>
	 * </ol>
	 */
	private boolean shouldOverrideFaceAssignment(Edge e, Edge current_edge,
			Vertex current_source, Vertex current_target, Vertex out_target,
			Pair<Vertex, Vertex> segment, Pair<Face, Face> nulled_faces, Pair<Face, Face> nullFace) {
		// Case 1: edge belongs to a face that has been nulled
		if (e.face.equals(nulled_faces.getFirst()) || e.face.equals(nulled_faces.getSecond())) {
			return true;
		}
		// Case 2: null-edge transition (only one null-edge in succession)
		if (e.type != EdgeType.NULLEDGE || current_edge.type == EdgeType.NULLEDGE) {
			return false;
		}
		// Must not be on a null-face
		if (e.face.equals(nullFace.getFirst()) || e.face.equals(nullFace.getSecond())) {
			return false;
		}
		// Valid transitions: sep→end, end→end, or target is a segment endpoint
		var sepToEnd = current_target.type == VertexType.SEPPOINT && out_target.type == VertexType.ENDPOINT;
		var endToEnd = current_source.type == VertexType.ENDPOINT && current_target.type == VertexType.ENDPOINT;
		var targetIsSegEndpoint = out_target.equals(segment.getFirst()) || out_target.equals(segment.getSecond());
		return sepToEnd || endToEnd || targetIsSegEndpoint;
	}

	/**
	 * Removes the IN vertices of the delete-tree; removes the IN vertices stored in
	 * {@code inVertices} (and associated IN-NEW edges)
	 */
	private void removeVertexSet() {
		for (Vertex v : inVertices) {
			// it should now be safe to delete all IN vertices
			assert (v.status == VertexStatus.IN) : "v.status == VertexStatus.IN";
			g.deleteVertex(v); // this also removes edges connecting to v
			modifiedVertices.remove(v);
		}
	}

	// remove all SPLIT type vertices on the HEFace \a f
	private void removeSplitVertex(Face f) {
		assert (checkFace(f)) : "face check failed";

		Vertex v;
		while ((v = findSplitVertex(f)) != null) {
			assert (v.type == VertexType.SPLIT) : "vertex must be SPLIT type for removal";

			g.removeDeg2Vertex(v);
			modifiedVertices.remove(v);

			assert (checkFace(f)) : "face check failed";
		}

		assert (checkFace(f)) : "face check failed";
	}

	// reset status of modified_vertices and incident_faces
	/**
	 * At the end of an incremental insertion of a new site, reset the status of
	 * modified_vertices to UNDECIDED and incident_faces to NONINCIDENT so that we
	 * are ready for the next insertion.
	 */
	private void resetStatus() {
		for (Vertex v : modifiedVertices) {
			v.resetStatus();
		}
		modifiedVertices.clear();
		for (Face f : incidentFaces) {
			f.setStatus(FaceStatus.NONINCIDENT);
		}
		incidentFaces.clear();
		inVertices.clear();
	}

	// count number of NEW vertices on the given face \a f
	private int numNewVertices(Face f) {
		var current = f.getEdge();
		var start = current;
		var count = 0;
		do {
			var v = current.target;
			if ((v.status == VertexStatus.NEW) && (v.type != VertexType.SEPPOINT)) {
				count++;
			}
			current = current.next;
		} while (!current.equals(start));
		return count;
	}
}
