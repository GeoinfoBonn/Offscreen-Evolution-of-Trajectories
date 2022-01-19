package de.geoinfoBonn.offscreenEvolution.evolution;

import java.awt.geom.Point2D;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;

import de.geoinfoBonn.offscreenEvolution.main.LookBackType;
import de.geoinfoBonn.offscreenEvolution.main.Variant;
import de.geoinfoBonn.offscreenEvolution.shapefinder.CakeslicePreparer;
import de.geoinfoBonn.offscreenEvolution.shapefinder.MinimalKShapeFinder;
import de.geoinfoBonn.offscreenEvolution.shapefinder.ShapePreparer;
import de.geoinfoBonn.offscreenEvolution.shapefinder.ShapePreparer.SignatureType;
import de.geoinfoBonn.offscreenEvolution.tracks.Track;
import de.geoinfoBonn.offscreenEvolution.tracks.TrackPoint;
import de.geoinfoBonn.offscreenEvolution.utils.GeometricTransformations;

public class ComputeEvolution {

	private static final Logger LOGGER = Logger.getLogger(ComputeEvolution.class.getName());
	private static final GeometryFactory GF = JTSFactoryFinder.getGeometryFactory();

	public static enum Direction {
		FORWARD, BACKWARD
	}

	public static LinkedList<EvolutionSymbol> computeEvolution(Direction direction, Track track, double[][] times,
			double kappa, TrackPoint current, TrackPoint previous, SignatureType type, Variant variant,
			double scaleFactor, LookBackType lbtype) {
		LinkedList<EvolutionSymbol> symbols = new LinkedList<>();
		int numShapes = times.length;

		// get shapes for each time span and store it in the list
		for (int iCake = 0; iCake < numShapes; iCake++) {
			List<TrackPoint> relevantPoints = getRelevantPoints(direction, track, times[iCake], current, variant,
					lbtype);
			// convert them to array
			if (relevantPoints == null || relevantPoints.isEmpty()) {
				LOGGER.severe("no relevant points found - choose a higher length");
				return null;
			}

			// Convert TrackPoints to Point2D
			ArrayList<Point2D> globalPoints = new ArrayList<>();
			for (TrackPoint p : relevantPoints)
				globalPoints.add(new Point2D.Double(p.getCoordinate().x, p.getCoordinate().y));

			// convert track points to Point2D
			Point2D currentPoint = new Point2D.Double(current.getCoordinate().x, current.getCoordinate().y);
			Point2D previousPoint = new Point2D.Double(previous.getCoordinate().x, previous.getCoordinate().y);

			// get parameters
			int n = globalPoints.size();
			int k = (int) Math.ceil(n * kappa);
			LOGGER.info("n = " + n + ", k = " + k);

			if (n == 1) {
				LOGGER.severe("only n=1 point found - very small cakeslice");
				return null;
			}

			// sort n points along ascending opening angle from last segment to point
			Collections.sort(globalPoints, new RadialPointComparator(currentPoint, previousPoint));

			// [smallest perimeter | index of corresponding point defining the right edge]
			double[] bestResult = new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
			// [r | phi] of optimal cakeslices
			double[] optParameters = new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };

			// iterate through list
			for (int l = 0; l <= globalPoints.size() - k; l++) {
				// current gamma defining right cakeslice edge
				double gamma = GeometricTransformations.angleFromX(currentPoint, globalPoints.get(l));

				// transform GPS points from global to local system
				ArrayList<Point2D> sites_loc = new ArrayList<>(globalPoints.size() - l);
				for (int j = l; j < globalPoints.size(); j++) {
					sites_loc.add(GeometricTransformations.global2local(gamma, currentPoint.getX(), currentPoint.getY(),
							globalPoints.get(j)));
				}

				// cakeslice
				switch (type) {
				case CAKE_SLICE:
					ShapePreparer cs = new CakeslicePreparer(true);
					double[] result_cake = MinimalKShapeFinder.findMinimalKShape(sites_loc, k, cs);

					// update the best result for perimeter if better solution is found
					double perimeter = cs.evaluateTargetFunction(result_cake[0], result_cake[1]);
					if (perimeter <= bestResult[0]) {
						bestResult = new double[] { perimeter, l };
						optParameters = result_cake;
					}
					break;
				default:
					LOGGER.severe("Type " + type + " is not implemented yet!");
					break;
				}
			}

			// return optimal right edge and optimal cakeslice
			LOGGER.fine("optimal right edge through point " + globalPoints.get((int) bestResult[1])
					+ " and optimal cakeslice parameter r = " + optParameters[0] + ", phi = "
					+ optParameters[1] * 180 / Math.PI + "[°]");

			EvolutionSymbol bestSymbol = new EvolutionSymbol.CakesliceSymbol(//
					new Coordinate(currentPoint.getX(), currentPoint.getY()), //
					optParameters[1], //
					optParameters[0], //
					GeometricTransformations.angleFromX(currentPoint, globalPoints.get((int) bestResult[1])), //
					scaleFactor, //
					relevantPoints.stream().map(p -> GF.createPoint(p.getCoordinate())).collect(Collectors.toList()));
			symbols.add(bestSymbol);
		}

		return symbols;
	}

	/**
	 * Depending on the direction and the selected variant, this function extracts
	 * all points from the track that are relevant for the creation for this shape.
	 * This means the next 'n' points are extracted (or in this case the points
	 * within the given time frame).
	 * 
	 * @param direction direction of the shape to calculate in relation to track
	 * @param track     input trajectory
	 * @param times     relevant time frame for the shape
	 * @param current   point where shape should be attached to
	 * @param variant   all or just intermediate points?
	 * @return list of all relevant track points
	 */
	private static List<TrackPoint> getRelevantPoints(Direction direction, Track track, double[] times,
			TrackPoint current, Variant variant, LookBackType type) {
		List<TrackPoint> n_points = null;

		if (direction == Direction.FORWARD) {
			switch (variant) {
			case ALL_POINTS:
				if (type == LookBackType.TIME)
					n_points = track.getPointsAfter(current.getTimestamp(), Duration.ofSeconds((int) times[1]));
				else if (type == LookBackType.DISTANCE)
					n_points = track.getPointsBetween(current, times[1]);
				else if (type == LookBackType.RANGE)
					n_points = track.getPointsInRange(current, times[1]);
				break;

			case INTERMEDIATE_POINTS:
				if (type == LookBackType.TIME) {
					// start of time interval (start at current)
					List<TrackPoint> startpoints = track.getPointsAfter(current.getTimestamp(),
							Duration.ofSeconds((int) times[0]));
					// last point of first point sequence = start point for actual interval
					TrackPoint start = current;
					if (!startpoints.isEmpty())
						start = startpoints.get(startpoints.size() - 1);

					// add to start point the length of the current time interval
					n_points = track.getPointsAfter(start.getTimestamp(),
							Duration.ofSeconds((int) (times[1] - times[0])));
				} else if (type == LookBackType.DISTANCE) {
					// start of time interval (start at current)
					List<TrackPoint> startpoints = track.getPointsBetween(current, times[0]);

					// last point of first point sequence = start point for actual interval
					TrackPoint start = current;
					if (!startpoints.isEmpty())
						start = startpoints.get(startpoints.size() - 1);

					// add to start point the length of the current time interval
					n_points = track.getPointsBetween(start, times[1] - times[0]);
				} else if (type == LookBackType.RANGE) {
					// start of time interval (start at current)
					List<TrackPoint> startpoints = track.getPointsInRange(current, times[0]);

					// last point of first point sequence = start point for actual interval
					TrackPoint start = current;
					if (!startpoints.isEmpty())
						start = startpoints.get(startpoints.size() - 1);

					// add to start point the length of the current time interval
					n_points = track.getPointsInRange(start, times[1] - times[0]);
				}
				break;

			case ONE_SIGNATURE:
				LOGGER.warning("DEPRECATED");
				break;

			default:
				break;
			}
		} else if (direction == Direction.BACKWARD) {
			switch (variant) {
			case ALL_POINTS:
				if (type == LookBackType.TIME)
					n_points = track.getPointsBefore(current.getTimestamp(), Duration.ofSeconds((int) times[1]));
				else if (type == LookBackType.DISTANCE)
					n_points = track.getPointsBetween(current, -times[1]);
				else if (type == LookBackType.RANGE)
					n_points = track.getPointsInRange(current, -times[1]);
				break;

			case INTERMEDIATE_POINTS:
				if (type == LookBackType.TIME) {
					List<TrackPoint> startpoints = track.getPointsBefore(current.getTimestamp(),
							Duration.ofSeconds((int) times[0]));
					// first point of first point sequence = start point for actual interval
					TrackPoint start = current;
					if (!startpoints.isEmpty())
						start = startpoints.get(0);

					// add to start point the length of the current time interval
					n_points = track.getPointsBefore(start.getTimestamp(),
							Duration.ofSeconds((int) (times[1] - times[0])));
				} else if (type == LookBackType.DISTANCE) {
					List<TrackPoint> startpoints = track.getPointsBetween(current, -times[0]);
					// first point of first point sequence = start point for actual interval
					TrackPoint start = current;
					if (!startpoints.isEmpty())
						start = startpoints.get(0);

					// add to start point the length of the current time interval
					n_points = track.getPointsBetween(start, -(times[1] - times[0]));
				} else if (type == LookBackType.RANGE) {
					List<TrackPoint> startpoints = track.getPointsInRange(current, -times[0]);
					// first point of first point sequence = start point for actual interval
					TrackPoint start = current;
					if (!startpoints.isEmpty())
						start = startpoints.get(0);

					// add to start point the length of the current time interval
					n_points = track.getPointsInRange(start, -(times[1] - times[0]));
				}
				break;

			case ONE_SIGNATURE:
				LOGGER.warning("DEPRECATED");
				break;

			default:
				break;
			}
		}

		return n_points;
	}
}
