package de.geoinfoBonn.offscreenEvolution.shapefinder;

import java.awt.geom.Point2D;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.util.GeometricShapeFactory;

import de.geoinfoBonn.offscreenEvolution.utils.GeometricTransformations;

public class CakeslicePreparer implements ShapePreparer {

	private boolean usePerimeter;

	public CakeslicePreparer() {
		this.usePerimeter = true;
	}

	public CakeslicePreparer(boolean usePerimeter) {
		this.usePerimeter = usePerimeter;
	}

	@Override
	public Point2D[] transformPoints(Point2D[] input) {
		// cartesian to polar
		Point2D[] polar = new Point2D[input.length];
		for (int i = 0; i < input.length; i++) {
			polar[i] = GeometricTransformations.cart2pol(input[i]);
		}
		return polar;
	}

	@Override
	public double evaluateTargetFunction(double r, double phi) {
		if (usePerimeter) {
//			 if (phi >= Math.PI / 2) {
//			 System.err.println("Phi zu gro�!");
//			 return Double.MAX_VALUE;
//			 }
			return 2 * r + r * phi;
		} else {
			// area
			return 0.5 * phi * Math.pow(r, 2);
		}
	}

	@Override
	public Polygon createShape(double r, double phi) {
		GeometryFactory gf = new GeometryFactory();
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory(gf);
		shapeFactory.setCentre(new Coordinate(0, 0));
		shapeFactory.setSize(2 * r);
		LineString arc = shapeFactory.createArc(0, phi);
		Coordinate[] coords = new Coordinate[arc.getCoordinates().length + 2];
		coords[0] = new Coordinate(0, 0);
		for (int i = 0; i < arc.getCoordinates().length; ++i) {
			coords[i + 1] = arc.getCoordinateN(i);
		}
		coords[coords.length - 1] = new Coordinate(0, 0);
		Polygon poly = gf.createPolygon(coords);
		return poly;
	}

}
