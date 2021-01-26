package org.rogach.jopenvoronoi;

/**
 * Denotes the permanent type of a vertex in the diagram.
 */
public enum VertexType {

	/**
	 * OUTER vertices are special vertices added in init(), should have degree==4
	 */
	OUTER,

	/**
	 * NORMAL are normal voronoi-vertices, should have degree==6 (degree 3 graph
	 * with double-edges)
	 */
	NORMAL,

	/**
	 * POINTSITE are point sites, should have degree==0
	 */
	POINTSITE,

	/**
	 * ENDPOINT vertices are end-points of line-segments or arc-segments
	 */
	ENDPOINT,

	/**
	 * separator start-vertices on a null-face
	 */
	SEPPOINT,

	/**
	 * APEX vertices split quadratic edges at their apex(closest point to site)
	 */
	APEX,

	/**
	 * split-vertices of degree==2 to avoid loops in the delete-tree
	 */
	SPLIT
};
