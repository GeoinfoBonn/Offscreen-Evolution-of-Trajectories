package de.geoinfoBonn.offscreenEvolution.utils;

import java.awt.geom.Point2D;

import org.locationtech.jts.geom.Coordinate;

public abstract class GeometricTransformations {

	/**
	 * 
	 * @param p local cartesian coordinate [x, y]
	 * @return polar coordinate [r, phi]
	 */
	public static Point2D cart2pol(Point2D p) {
		double r = Math.sqrt(p.getX() * p.getX() + p.getY() * p.getY());
		double phi = java.lang.Math.atan2(p.getY(), p.getX());
		Point2D polar = new Point2D.Double(r, phi);
		return polar;
	}

	/**
	 * 
	 * @param p local polar coordinate [r, phi]
	 * @return cartesian coordinate [x, y]
	 */
	public static Point2D pol2cart(Point2D p) {
		double r = p.getX();
		double phi = p.getY();
		double x = r * Math.cos(phi);
		double y = r * Math.sin(phi);
		Point2D cart = new Point2D.Double(x, y);
		return cart;
	}

	/**
	 * 
	 * @param xi  polar coordinates [xi, phi] of triangle
	 * @param phi polar coordinates [xi, phi] of triangle [rad]
	 * @return cartesian coordinates [x, y] of triangle
	 */
	public static Point2D trianglePol2cart(double xi, double phi) {
		double x = xi;
		double y = xi * Math.tan(phi);
		return new Point2D.Double(x, y);
	}

	/**
	 * 
	 * @param origin_glob point given in global [x, y] and defines local xi-axis
	 * @param end_glob    point given in global [x, y] and defines local xi-axis
	 * @param p_glob      point to be transformed; given in global [x, y]
	 * @return local point [xi, eta] relative to local xi-axis defined by
	 *         origin_glob and end_glob
	 */
	public static Point2D global2local(Point2D origin_glob, Point2D end_glob, Point2D p_glob) {
		double gamma = angleFromX(origin_glob, end_glob);
		double tx = origin_glob.getX(); // point in which triangle has to be created
		double ty = origin_glob.getY();
		double xi = Math.cos(-gamma) * (p_glob.getX() - tx) - Math.sin(-gamma) * (p_glob.getY() - ty);
		double eta = Math.sin(-gamma) * (p_glob.getX() - tx) + Math.cos(-gamma) * (p_glob.getY() - ty);
		return new Point2D.Double(xi, eta);
	}

	/**
	 * 
	 * @param gamma  angle of orientation between global x-axis and local xi-axis;
	 *               orientation counterclockwise
	 * @param tx     global x-coordinate of origin point; translation in x-direction
	 * @param ty     global y-coordinate of origin point; translation in y-direction
	 * @param p_glob global point [x, y] which has to be transformed
	 * @return local point [xi, eta]
	 */
	public static Point2D global2local(double gamma, double tx, double ty, Point2D p_glob) {
		double xi = Math.cos(-gamma) * (p_glob.getX() - tx) - Math.sin(-gamma) * (p_glob.getY() - ty);
		double eta = Math.sin(-gamma) * (p_glob.getX() - tx) + Math.cos(-gamma) * (p_glob.getY() - ty);
		return new Point2D.Double(xi, eta);
	}

	/**
	 * 
	 * @param gamma known parameter: angle of orientation between global x-axis and
	 *              local xi-axis; orientation counterclockwise
	 * @param tx    known parameter: global x-coordinate of origin point;
	 *              translation in x-direction
	 * @param ty    known parameter: global y-coordinate of origin point;
	 *              translation in y-direction
	 * @param p_loc local point [xi, eta] which has to be transformed
	 * @return global point [x, y]
	 */
	public static Point2D local2global(double gamma, double tx, double ty, Point2D p_loc) {
		double x = tx + Math.cos(gamma) * p_loc.getX() - Math.sin(gamma) * p_loc.getY();
		double y = ty + Math.sin(gamma) * p_loc.getX() + Math.cos(gamma) * p_loc.getY();
		return new Point2D.Double(x, y);
	}

	/**
	 * 
	 * @param p_loc local polar coordinate [r, phi] (respectively [xi, phi] for
	 *              triangle)
	 * @return false if point lies behind (means: phi > 90 or phi < -90)
	 */
	public static boolean checkPointLiesBehind(Point2D p_loc) {
		double phi = p_loc.getY();
		if (phi > Math.PI / 2 || phi < -Math.PI / 2) {
			return false; // point has to be sorted out
		} else {
			return true; // point can be used
		}
	}

	/**
	 * 
	 * @param p_loc local polar coordinate [r, phi] (respectively [xi, phi] for
	 *              triangle)
	 * @return local polar point with positive phi; lies left-hand side of the local
	 *         xi-axis
	 */
	public static Point2D flipPointAlongXiAxis(Point2D p_loc) {
		double phi = p_loc.getY();
		if (phi < 0) {
			double newPhi = Math.abs(phi);
			return new Point2D.Double(p_loc.getX(), newPhi);
		} else {
			return new Point2D.Double(p_loc.getX(), p_loc.getY());
		}
	}

	/**
	 * 
	 * @param p Point2D, point where the new point q is going to be attached
	 * @param r distance to new point q [m]
	 * @param t Richtungswinkel to new point q [rad]
	 * @return q polarly attached new point
	 */
	public static Point2D polarAttachment(Point2D p, double r, double t) {
		double px = p.getX();
		double py = p.getY();
		double qx = px + r * Math.sin(t);
		double qy = py + r * Math.cos(t);
		Point2D q = new Point2D.Double(qx, qy);
		return q;
	}

	/**
	 * 
	 * @param gamma positive orientation angle from horizontal axis
	 *              counter-clockwise
	 * @return Richtungswinkel t from north-orientated vertical axis clockwise
	 */
	public static double richtungswinkelGivenGamma(double gamma) {
		if (gamma <= Math.PI / 2) {
			return Math.PI / 2 - gamma;
		}
		if (gamma <= Math.PI) {
			return (2 * Math.PI) - (gamma - Math.PI / 2);
		}
		if ((gamma >= Math.PI) && (gamma <= 3 * Math.PI / 2)) {
			return (2 * Math.PI) - (gamma - Math.PI / 2);
		}
		if (gamma <= 2 * Math.PI) {
			return Math.PI / 2 + (2 * Math.PI - gamma);
		}
		if (gamma > 2 * Math.PI) {
			return (2 * Math.PI - gamma) + Math.PI / 2;
		}
		return 0.0;

	}

	/**
	 * Calculates the inner angle spanned by the two line segments
	 * anchorPoint->referencePoint and anchorPoint->spanningPoint. The returned
	 * angle is given counter-clockwise in radian and normalized to the interval
	 * [0,2pi).
	 * 
	 * @param referencePoint the angle is measured beginning at this point
	 * @param anchorPoint    the angle is measured around this point
	 * @param spanningPoint  the angle is spanned by this point
	 * @return
	 */
	public static double innerAngle(Point2D referencePoint, Point2D anchorPoint, Point2D spanningPoint) {
		if (anchorPoint.equals(referencePoint) || anchorPoint.equals(spanningPoint))
			throw new IllegalArgumentException(
					"inner angle not defined when one of the points is equal to the anchor!");

		double gammaR = angleFromX(anchorPoint, referencePoint);
		double gammaS = angleFromX(anchorPoint, spanningPoint);

		double phi = gammaS - gammaR;
		if (phi < 0)
			phi += 2 * Math.PI;
		return phi;
	}

	/**
	 * Returns the angle between the x-axis and the line segment
	 * anchorPoint->spanningPoint. The returned angle is given in radian and
	 * normalized to [0;2pi).
	 * 
	 * @param anchorPoint   the angle is measured around this point
	 * @param spanningPoint the angle is spanned from the x-axis to this point
	 * @return angle gamma [rad] between global x-axis and line through
	 *         {@link anchorPoint} and {@link spanningPoint}; orientation
	 *         counterclockwise
	 */
	public static double angleFromX(Point2D anchorPoint, Point2D spanningPoint) {
		if (anchorPoint.equals(spanningPoint))
			throw new IllegalArgumentException("angle between two identical points is not defined!");

		double gamma = Math.atan2(spanningPoint.getY() - anchorPoint.getY(), spanningPoint.getX() - anchorPoint.getX());
		if (gamma < 0)
			gamma += 2 * Math.PI;
		return gamma;
	}

	public static double angleFromX(Coordinate anchorPoint, Coordinate spanningPoint) {
		if (anchorPoint.equals(spanningPoint))
			throw new IllegalArgumentException("angle between two identical points is not defined!");

		double gamma = Math.atan2(spanningPoint.getY() - anchorPoint.getY(), spanningPoint.getX() - anchorPoint.getX());
		if (gamma < 0)
			gamma += 2 * Math.PI;
		return gamma;
	}
}
