package org.rogach.jopenvoronoi.offset;

import java.util.ArrayList;
import java.util.List;

//a single offset loop
public class OffsetLoop {
	/** list of offsetvertices in this loop */
	public List<OffsetVertex> vertices = new ArrayList<>();
	public double offsetDistance;

	void add(OffsetVertex v) {
		vertices.add(v);
	}
}
