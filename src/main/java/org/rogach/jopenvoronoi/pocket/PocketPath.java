package org.rogach.jopenvoronoi.pocket;

import java.util.ArrayList;
import java.util.List;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * Converts a sequence of {@link MIC Maximal Inscribed Circles} into a
 * continuous toolpath of line and arc segments.
 * <p>
 * For each consecutive pair of MICs the path consists of:
 * <ol>
 *   <li>An <b>arc</b> on the previously-cut circle connecting the two
 *       bi-tangent contact points on that circle (tracing the "floor" that
 *       is already machined).</li>
 *   <li>A <b>line</b> along one external bi-tangent from the old circle to
 *       the new circle.</li>
 *   <li>An <b>arc</b> on the new circle connecting its two bi-tangent
 *       contact points (cutting fresh material).</li>
 *   <li>A <b>line</b> along the other external bi-tangent back toward the
 *       old circle.</li>
 * </ol>
 * The result is a zig-zag / spiral-like path that sweeps material between
 * successive MICs. When a {@linkplain MIC#newBranch new branch} is
 * encountered, a rapid/retract move is inserted to reposition the tool.
 *
 * <p><b>Usage</b></p>
 * <pre>{@code
 * MedialAxisPocket pocket = new MedialAxisPocket(vd.getDiagram());
 * pocket.setWidth(0.05);
 * pocket.run();
 *
 * for (List<MIC> comp : pocket.getMicComponents()) {
 *     List<PocketPath.Segment> path = PocketPath.toPath(comp);
 *     for (PocketPath.Segment seg : path) {
 *         if (seg.isArc()) {
 *             // arc move to seg.end around seg.center, radius seg.radius
 *         } else {
 *             // line move to seg.end
 *         }
 *     }
 * }
 * }</pre>
 *
 * @see MIC
 * @see MedialAxisPocket
 */
public final class PocketPath {

	private PocketPath() {
	}

	/**
	 * A single segment of the pocketing toolpath — either a straight line or a
	 * circular arc.
	 */
	public static class Segment {

		/** start point of this segment */
		public final Point start;
		/** end point of this segment */
		public final Point end;
		/**
		 * center of the arc, or {@code null} for a line segment
		 */
		public final Point center;
		/**
		 * radius of the arc (positive), or {@code -1} for a line segment
		 */
		public final double radius;
		/** {@code true} for clockwise arcs */
		public final boolean clockwise;
		/**
		 * Segment type.
		 *
		 * @see Type
		 */
		public final Type type;

		private Segment(Point start, Point end, Point center, double radius, boolean clockwise, Type type) {
			this.start = start;
			this.end = end;
			this.center = center;
			this.radius = radius;
			this.clockwise = clockwise;
			this.type = type;
		}

		/** Returns {@code true} if this segment is an arc. */
		public boolean isArc() {
			return type == Type.ARC;
		}

		/** Returns {@code true} if this segment is a straight line. */
		public boolean isLine() {
			return type == Type.LINE;
		}

		/** Returns {@code true} if this segment is a rapid/retract move. */
		public boolean isRapid() {
			return type == Type.RAPID;
		}

		@Override
		public String toString() {
			if (isArc()) {
				return String.format("Arc(%s → %s, c=%s, r=%.4f, %s)", start, end, center, radius,
						clockwise ? "CW" : "CCW");
			} else if (isRapid()) {
				return String.format("Rapid(%s → %s)", start, end);
			} else {
				return String.format("Line(%s → %s)", start, end);
			}
		}
	}

	/**
	 * Segment type.
	 */
	public enum Type {
		/** A straight cutting move. */
		LINE,
		/** A circular arc cutting move. */
		ARC,
		/**
		 * A rapid (non-cutting) repositioning move, emitted when the tool
		 * must jump to a new medial-axis branch.
		 */
		RAPID
	}

	/**
	 * Converts a list of MICs (one connected component) into a continuous
	 * toolpath.
	 * <p>
	 * The first MIC in the list is the starting circle. For each subsequent
	 * MIC the method emits:
	 * <ol>
	 *   <li>If the MIC is a {@linkplain MIC#newBranch new branch}, a
	 *       {@link Type#RAPID RAPID} segment from the current position to the
	 *       branch origin.</li>
	 *   <li>An arc on the previously-cut circle between bi-tangent points
	 *       {@link MIC#t1} and {@link MIC#t2}.</li>
	 *   <li>A line along the bi-tangent from {@link MIC#t1} to
	 *       {@link MIC#t3}.</li>
	 *   <li>An arc on the new circle between bi-tangent points
	 *       {@link MIC#t3} and {@link MIC#t4}.</li>
	 *   <li>A line along the bi-tangent from {@link MIC#t4} back to
	 *       {@link MIC#t2}.</li>
	 * </ol>
	 *
	 * @param mics MIC list from one connected component (as returned by
	 *             {@link MedialAxisPocket#getMicComponents()})
	 * @return ordered list of toolpath segments
	 */
	public static List<Segment> toPath(List<MIC> mics) {
		List<Segment> path = new ArrayList<>();
		if (mics == null || mics.isEmpty()) {
			return path;
		}

		Point cursor = mics.get(0).c2; // start at center of first (largest) MIC

		for (int i = 1; i < mics.size(); i++) {
			MIC mic = mics.get(i);
			if (mic.t1 == null || mic.t2 == null || mic.t3 == null || mic.t4 == null) {
				// degenerate MIC (c1 == c2), skip
				continue;
			}

			if (mic.newBranch) {
				// rapid to the branch origin
				path.add(rapid(cursor, mic.cPrev));
				cursor = mic.cPrev;
			}

			// 1) arc on previous circle: from t1 to t2 around c1
			boolean cw1 = isClockwise(mic.t1, mic.c1, mic.t2);
			path.add(arc(mic.t1, mic.t2, mic.c1, mic.r1, cw1));

			// 2) line along bi-tangent: t1 → t3
			path.add(line(mic.t1, mic.t3));

			// 3) arc on new circle: from t3 to t4 around c2
			boolean cw2 = isClockwise(mic.t3, mic.c2, mic.t4);
			path.add(arc(mic.t3, mic.t4, mic.c2, mic.r2, cw2));

			// 4) line along bi-tangent: t4 → t2
			path.add(line(mic.t4, mic.t2));

			cursor = mic.t2;
		}
		return path;
	}

	/**
	 * Determines the winding direction of an arc from {@code start} to
	 * {@code end} around {@code center}.
	 */
	private static boolean isClockwise(Point start, Point center, Point end) {
		// cross product of (start−center) × (end−center)
		double ax = start.x - center.x;
		double ay = start.y - center.y;
		double bx = end.x - center.x;
		double by = end.y - center.y;
		return (ax * by - ay * bx) < 0;
	}

	private static Segment line(Point start, Point end) {
		return new Segment(start, end, null, -1, false, Type.LINE);
	}

	private static Segment arc(Point start, Point end, Point center, double radius, boolean cw) {
		return new Segment(start, end, center, radius, cw, Type.ARC);
	}

	private static Segment rapid(Point start, Point end) {
		return new Segment(start, end, null, -1, false, Type.RAPID);
	}
}
