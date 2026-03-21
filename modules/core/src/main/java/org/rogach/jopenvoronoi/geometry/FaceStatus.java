package org.rogach.jopenvoronoi.geometry;

public enum FaceStatus {

	/**
	 * INCIDENT faces contain one or more IN-vertex
	 */
	INCIDENT,
	/**
	 * NONINCIDENT faces contain only OUT/UNDECIDED-vertices
	 */
	NONINCIDENT
}
