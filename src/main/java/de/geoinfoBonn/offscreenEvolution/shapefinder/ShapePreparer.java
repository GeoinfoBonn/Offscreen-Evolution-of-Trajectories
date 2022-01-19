package de.geoinfoBonn.offscreenEvolution.shapefinder;

import java.awt.geom.Point2D;

import org.locationtech.jts.geom.Polygon;

public interface ShapePreparer {

	public static enum SignatureType {
		RECTANGLE, TRIANGLE, CAKE_SLICE
	}

	public Point2D[] transformPoints(Point2D[] input); // prepare local coordinates i.e. polar

	public double evaluateTargetFunction(double a, double b); // what has to be minimized i.e. perimeter

	public Polygon createShape(double a, double b); // for easy visualization and tests

}
