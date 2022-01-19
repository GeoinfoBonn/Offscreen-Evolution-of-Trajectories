package de.geoinfoBonn.offscreenEvolution.shapefinder;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

public class MinimalKShapeFinder {

	/**
	 * 
	 * @param localPoints local cartesian coordinates xi, eta
	 * @param k           amount of points in shape
	 * @param sp          sets target shape, prepares the input points
	 * @return local polar coordinates [xi, phi] of minimal shape
	 */
	public static double[] findMinimalKShape(Point2D[] localPoints, int k, ShapePreparer sp) {

		// prepare the local points, i.e. transform them to polar coordinates
		// (no mirroring or sorting because all points ie left to the right edge)
		Point2D[] points = sp.transformPoints(localPoints);

		if (points.length < k) {
			System.err.println("Too less sites around here to use!");
			return new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
		}

		// array sort in x- and y-direction
		Point2D[] X = MinimalKShapeFinder.sortX(points);
		Point2D[] Y = MinimalKShapeFinder.sortY(points);

		// index in array Y of point with highest y-coordinate among the k start-points
		int i_Y = MinimalKShapeFinder.indexHighestY(X, k, Y);

		// currently known optimal parameters
		double x_curr = X[k - 1].getX();
		double y_curr = Y[i_Y].getY();
		double x_best = x_curr;
		double y_best = y_curr;
		double opt = sp.evaluateTargetFunction(x_best, y_best);

		// iterations increasing x
		for (int i_X = k; i_X < points.length; i_X++) {
			// do not consider point if its y-value is greater than the current value
			if (X[i_X].getY() >= y_curr)
				continue;

			// update x-parameter of the rec and iterate over decreasing y
			x_curr = X[i_X].getX();

			// iterations decreasing y
			for (i_Y--; i_Y >= k; i_Y--) {
				// do not consider point if its x-value is greater than the current value
				if (Y[i_Y].getX() > x_curr)
					continue;

				// if yes: update h-parameter of opt. rec
				y_curr = Y[i_Y].getY();

				// current target result of opt. rec using interface
				double target = sp.evaluateTargetFunction(x_curr, y_curr);
				// update opt. target but avoid deformations (no zero)
				if (target <= opt && target != 0.0) {
					opt = target;
					y_best = y_curr;
					x_best = x_curr;
				}
				break;
			}
		}
		double[] result = new double[] { x_best, y_best };
		return result;
	}

	/**
	 * Sort points ascending in the x-coordinate (and descending in the y-coordinate
	 * in case of ties).
	 * 
	 * @param points input point array
	 * @return point array sorted by x-coordinate
	 */
	private static Point2D[] sortX(Point2D[] points) {
		Point2D[] X = new Point2D[points.length];
		for (int i = 0; i < points.length; i++)
			X[i] = points[i];

		// sort x ascending by x and descending by y
		Arrays.sort(X, new Comparator<Point2D>() {
			@Override
			public int compare(Point2D a, Point2D b) {
				int xComp = Double.compare(a.getX(), b.getX());
				if (xComp == 0) {
					return -Double.compare(a.getY(), b.getY());
				} else {
					return xComp;
				}
			}
		});

		return X;
	}

	/**
	 * Sort points ascending in the y-coordinate (and descending in the x-coordinate
	 * in case of ties).
	 * 
	 * @param points input point array
	 * @return point array sorted by y-coordinate
	 */
	private static Point2D[] sortY(Point2D[] points) {
		Point2D[] Y = new Point2D[points.length];
		for (int i = 0; i < points.length; i++)
			Y[i] = points[i];

		// sort y ascending by y and descending by x
		Arrays.sort(Y, new Comparator<Point2D>() {
			@Override
			public int compare(Point2D a, Point2D b) {
				int yComp = Double.compare(a.getY(), b.getY());
				if (yComp == 0) {
					return -Double.compare(a.getX(), b.getX());
				} else {
					return yComp;
				}
			}
		});

		return Y;
	}

	/**
	 * Searches the point with the largest y-coordinate in the first {@code k}
	 * values of {@code X} and returns its index in {@code Y}.
	 * 
	 * @param X
	 * @param Y
	 * @param k
	 * @return index in Y
	 */
	private static int indexHighestY(Point2D[] X, int k, Point2D[] Y) {
		// find index of highest point in array Y of the first k points in array X (->
		// not globally)
		int indexHighestY = 0;
		double aktMaxY = 0;
		for (int i = 0; i < k; ++i) {
			Point2D aktPoint = X[i];
			if (aktPoint.getY() > aktMaxY) {
				aktMaxY = aktPoint.getY();
				indexHighestY = i;
			}
		}
		// find corresponding index of indexHighestY (which is until here reffered to
		// array X) in array Y
		Point2D point = X[indexHighestY];
		for (int i = 0; i < Y.length; i++) {
			if (Y[i].equals(point)) {
				indexHighestY = i;
				break;
			}
		}
		return indexHighestY;
	}

	/**
	 * 
	 * @param localPoints local cartesian coordinates xi, eta
	 * @param k           amount of points in shape
	 * @param sp          sets target shape, prepares the input points
	 * @return local polar coordinates [xi, phi] of minimal shape
	 */
	public static <P extends Point2D> double[] findMinimalKShape(ArrayList<P> localPoints, int k, ShapePreparer sp) {
		// convert given List to normal array
		Point2D[] localArray = new Point2D[localPoints.size()];
		for (int i = 0; i < localPoints.size(); i++) {
			localArray[i] = new Point2D.Double(localPoints.get(i).getX(), localPoints.get(i).getY());
		}

		return MinimalKShapeFinder.findMinimalKShape(localArray, k, sp);
	}

	/**
	 * Small test instance. Visualization can be found in
	 * offscreenevolution-paper/docs/algorithm_animation.ipe
	 * 
	 * @author axel
	 * @param args
	 */
	public static void main(String[] args) {
		ArrayList<Point2D> points = new ArrayList<>();
		points.add(new Point2D.Double(21, 11));
		points.add(new Point2D.Double(10, 18));
		points.add(new Point2D.Double(24, 4));
		points.add(new Point2D.Double(7, 4));
		points.add(new Point2D.Double(25, 18));
		points.add(new Point2D.Double(4, 12));
		points.add(new Point2D.Double(25, 8));
		points.add(new Point2D.Double(3, 17));
		points.add(new Point2D.Double(14, 8));
		points.add(new Point2D.Double(2, 2));
		points.add(new Point2D.Double(19, 15));
		points.add(new Point2D.Double(19, 5));
		points.add(new Point2D.Double(11, 11));
		points.add(new Point2D.Double(4, 6));
		points.add(new Point2D.Double(14, 2));

		int k = 10;

		double[] res = MinimalKShapeFinder.findMinimalKShape(points, k, new RectanglePreparer());
		System.out.println("Result (perimeter): x=" + res[0] + ", y=" + res[1]);

		res = MinimalKShapeFinder.findMinimalKShape(points, k, new RectanglePreparer(false));
		System.out.println("Result (area):      x=" + res[0] + ", y=" + res[1]);
	}
}
