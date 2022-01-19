package de.geoinfoBonn.offscreenEvolution.shapefinder;

import java.awt.geom.Point2D;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;

public class RectanglePreparer implements ShapePreparer {

	private boolean usePerimeter;

	public RectanglePreparer() {
		this.usePerimeter = true;
	}

	public RectanglePreparer(boolean usePerimeter) {
		this.usePerimeter = usePerimeter;
	}

	@Override
	public Point2D[] transformPoints(Point2D[] input) {
		Point2D[] ret = new Point2D[input.length];
		for (int i = 0; i < input.length; i++) {
			ret[i] = new Point2D.Double(input[i].getX(), input[i].getY());
		}
		return ret;
	}

	@Override
	public double evaluateTargetFunction(double h, double w) {
		if (usePerimeter)
			return 2 * h + 2 * w;
		return h * w;
	}

	@Override
	public Polygon createShape(double h, double w) {
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		Coordinate[] coords = new Coordinate[] { new Coordinate(0, 0), new Coordinate(h, 0), new Coordinate(h, w),
				new Coordinate(0, w), new Coordinate(0, 0) };
		Polygon poly = gf.createPolygon(coords);
		return poly;
	}

}
