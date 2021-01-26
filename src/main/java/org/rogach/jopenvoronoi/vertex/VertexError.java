package org.rogach.jopenvoronoi.vertex;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.rogach.jopenvoronoi.HalfEdgeDiagram;
import org.rogach.jopenvoronoi.geometry.Edge;
import org.rogach.jopenvoronoi.geometry.EdgeType;
import org.rogach.jopenvoronoi.geometry.Point;
import org.rogach.jopenvoronoi.site.Site;

//\brief error functor for edge-based desperate solver
///
//minimize error by searching for a point on the solution-edge
public class VertexError implements UnivariateFunction {
	HalfEdgeDiagram g; // < vd-graph
	Edge edge; // < existing edge on which we have positioned a new vertex
	Site s3; // < newly inserted Site

	// \param gi vd-graph
	// \param sln_edge solution edge
	// \param si3 newly inserted Site
	public VertexError(HalfEdgeDiagram gi, Edge sln_edge, Site si3) {
		this.g = gi;
		this.edge = sln_edge;
		this.s3 = si3;
	}

	// return the vertex-error t-d3 where
	// t3 is the distance from edge-point(t) to s3, and
	// t is the offset-distance of the solution
	@Override
	public double value(double t) {
		var p = edge_point(t);
		var s3_dist = p.sub(s3.apex_point(p)).norm();
		return Math.abs(t - s3_dist);
	}

	// return a point on the edge at given offset-distance
	// \param t offset-distance ( >= 0 )
	Point edge_point(double t) {
		Point p;
		if (edge.type == EdgeType.LINELINE) { // this is a workaround because the LINELINE edge-parameters are wrong? at
												// least in some cases?
			var src = edge.source;
			var trg = edge.target;
			var src_p = src.position;
			var trg_p = trg.position;
			var src_t = src.dist();
			var trg_t = trg.dist();
			// edge is src_p -> trg_p
			if (trg_t > src_t) {
				var frac = (t - src_t) / (trg_t - src_t);
				p = src_p.add(trg_p.sub(src_p).mult(frac));
			} else {
				var frac = (t - trg_t) / (src_t - trg_t);
				p = trg_p.add(src_p.sub(trg_p).mult(frac));
			}

		} else {
			p = edge.point(t);
		}
		return p;
	}
}
