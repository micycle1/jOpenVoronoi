package org.rogach.jopenvoronoi.geometry;

import org.rogach.jopenvoronoi.site.Site;

public class Face {
	public Edge edge;
	public Site site;
	public FaceStatus status;
	public boolean is_null_face;

	public Face() {
	}

	@Override
	public String toString() {
		var sb = new StringBuilder();
		sb.append("F(");
		var current = edge;
		var c = 0;
		do {
			if (current == null) {
				break;
			}
			sb.append(current.source.position);
			sb.append(">");
			current = current.next;
			c++;
		} while (current != edge && c < 100);
		if (c >= 100) {
			sb.append("...");
		}
		sb.append(")");
		return sb.toString();
	}
}
