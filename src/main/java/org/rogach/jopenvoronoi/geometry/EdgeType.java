package org.rogach.jopenvoronoi.geometry;

public enum EdgeType {

	/**
	 * Line edge between PointSite and PointSite
	 */
	LINE,
	/**
	 * Line edge between LineSite and LineSite
	 */
	LINELINE,
	/**
	 * Line edge between LineSite and LineSite (parallel case)
	 */
	PARA_LINELINE,
	/**
	 * special outer edge set by initialize()
	 */
	OUTEDGE,
	/**
	 * Parabolic edge between PointSite and LineSite
	 */
	PARABOLA,
	/**
	 * Separator edge between PointSite (endpoint) and LineSite or ArcSite
	 */
	ELLIPSE, HYPERBOLA, SEPARATOR,
	/**
	 * zero-length null-edge around a PointSite which is and endpoint
	 */
	NULLEDGE,
	/**
	 * pseudo-edge corresponding to a LineSite
	 */
	LINESITE,
	/**
	 * pseudo-edge corresponding to a LineSite
	 */
	ARCSITE

}
