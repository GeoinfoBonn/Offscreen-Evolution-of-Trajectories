package de.geoinfoBonn.offscreenEvolution.evolution;

import java.awt.geom.Point2D;
import java.util.Comparator;

import de.geoinfoBonn.offscreenEvolution.utils.GeometricTransformations;

public class RadialPointComparator implements Comparator<Point2D> {

	private Point2D anchorPoint; // defines rotating point of angle
	private Point2D referencePoint; // defines direction of zero degree

	public RadialPointComparator(Point2D anchorPoint, Point2D referencePoint) {
		this.anchorPoint = anchorPoint;
		this.referencePoint = referencePoint;
	}

	@Override
	public int compare(Point2D a, Point2D b) {
		if (a.equals(anchorPoint)) {
			return 1;
		} else if (b.equals(anchorPoint)) {
			return -1;
		}

		double angle_a = GeometricTransformations.innerAngle(referencePoint, anchorPoint, a);
		double angle_b = GeometricTransformations.innerAngle(referencePoint, anchorPoint, b);
		int c = Double.compare(angle_a, angle_b);
		return c;
	}

}
