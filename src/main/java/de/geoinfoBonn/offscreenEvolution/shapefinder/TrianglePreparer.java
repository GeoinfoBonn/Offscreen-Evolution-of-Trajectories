package de.geoinfoBonn.offscreenEvolution.shapefinder;

import java.awt.geom.Point2D;
import java.util.ArrayList;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;

import de.geoinfoBonn.offscreenEvolution.utils.GeometricTransformations;

public class TrianglePreparer implements ShapePreparer {

	private boolean usePerimeter = true;

	// public TrianglePreparer() {}

	public TrianglePreparer(boolean usePerimeter) {
		this.usePerimeter = usePerimeter;
	}

	@Override
	public Point2D[] transformPoints(Point2D[] input) {
		// cartesian to polar
		Point2D[] polar = new Point2D[input.length];
		for (int i = 0; i < input.length; i++) {
			Point2D temp = GeometricTransformations.cart2pol(input[i]);
			// set r = xi, only take phi
			polar[i] = new Point2D.Double(input[i].getX(), temp.getY());
		}
		// return polar; // old :-)
		// test phi
		ArrayList<Point2D> chosenPoints = new ArrayList<Point2D>(polar.length);
		for (int i = 0; i < polar.length; i++) {
			if (GeometricTransformations.checkPointLiesBehind(polar[i])) {
				Point2D chosenPoint = new Point2D.Double(polar[i].getX(), polar[i].getY());
				chosenPoints.add(chosenPoint); // use point
				// System.out.println("Point " + i+ " can be used.");
			} else {
				// System.out.println("Point " + i+ " cannot be used.");
			}
		}
		// if needed: flip point
		Point2D[] allPoints = new Point2D[chosenPoints.size()];
		for (int i = 0; i < chosenPoints.size(); i++) {
			allPoints[i] = GeometricTransformations.flipPointAlongXiAxis(chosenPoints.get(i));
		}
		// System.out.println("Number of points that can be used: " + allPoints.length);
		return allPoints;
	}

	@Override
	public double evaluateTargetFunction(double x, double phi) {
		if (usePerimeter) {
			if (phi >= Math.PI / 2) {
				System.err.println("Phi too large!");
				return Double.MAX_VALUE;
			}
			double y = Math.tan(phi) * x;
			double r = Math.sqrt(x * x + y * y);
			return x + y + r;
		} else {
			// area
			return 0.5 * Math.pow(x, 2) * Math.tan(phi);
		}

	}

	@Override
	public Polygon createShape(double xi, double phi) {
		Point2D cartesian = GeometricTransformations.trianglePol2cart(xi, phi);
		double x_tri = cartesian.getX();
		double y_tri = cartesian.getY();
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		Coordinate[] coords = new Coordinate[] { new Coordinate(0, 0), new Coordinate(x_tri, 0),
				new Coordinate(x_tri, y_tri), new Coordinate(0, 0) };
		Polygon poly = gf.createPolygon(coords);
		return poly;
	}

}
