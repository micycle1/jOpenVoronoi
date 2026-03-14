package org.rogach.jopenvoronoi.offset;

import org.rogach.jopenvoronoi.geometry.Point;

/**
 * Base class for offset elements.
 * <p>
 * Preliminary offset presentations. Experimental.
 */
public abstract class Ofs {
	/**
	 * Returns the radius of this offset element.
	 *
	 * @return radius, or {@code -1} if this element is a line
	 */
	public abstract double radius(); // {return -1;}

	/**
	 * Returns the center point of this offset element.
	 *
	 * @return center point for an arc element
	 */
	public abstract Point center(); // {return Point(0,0);}

	/**
	 * Returns the start point of this offset element.
	 *
	 * @return start point
	 */
	public abstract Point start();// {return Point(0,0);}

	/**
	 * Returns the end point of this offset element.
	 *
	 * @return end point
	 */
	public abstract Point end(); // {return Point(0,0);}
};
