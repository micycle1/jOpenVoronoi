package org.rogach.jopenvoronoi.offset;

import java.util.ArrayList;
import java.util.List;

//a single offset loop
public class OffsetLoop {
	public List<OffsetVertex> vertices = new ArrayList<>(); // < list of offsetvertices in this loop
	public double offset_distance;

	public void add(OffsetVertex v) {
		vertices.add(v);
	}
}
