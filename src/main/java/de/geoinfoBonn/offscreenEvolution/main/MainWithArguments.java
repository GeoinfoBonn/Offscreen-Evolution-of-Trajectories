package de.geoinfoBonn.offscreenEvolution.main;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

import de.geoinfoBonn.offscreenEvolution.evolution.ComputeEvolution;
import de.geoinfoBonn.offscreenEvolution.evolution.ComputeEvolution.Direction;
import de.geoinfoBonn.offscreenEvolution.evolution.EvolutionSymbol;
import de.geoinfoBonn.offscreenEvolution.io.Feature;
import de.geoinfoBonn.offscreenEvolution.io.FeatureReader;
import de.geoinfoBonn.offscreenEvolution.shapefinder.ShapePreparer.SignatureType;
import de.geoinfoBonn.offscreenEvolution.tracks.Track;
import de.geoinfoBonn.offscreenEvolution.tracks.TrackManager;
import de.geoinfoBonn.offscreenEvolution.tracks.TrackPoint;
import de.geoinfoBonn.offscreenEvolution.tracks.TrackPoint.PointType;
import de.geoinfoBonn.offscreenEvolution.utils.Tools;

public class MainWithArguments {

	private static final Logger LOGGER = Logger.getLogger(MainWithArguments.class.getName());

	public static void main(String args[]) throws IOException {
		LOGGER.info("Started offscreen-evaluation app.");

		RunConfig config = RunConfig.parseArguments(args);
		config.logRunConfig();

		if (!config.getOutputDir().exists())
			config.getOutputDir().mkdirs();

		// ################## read the trajectories from shapefile ##################
		List<Track> tracks = Track.importFromShapefile(config.getTrackFile().getAbsolutePath(), config.getIdColumn(),
				config.getSegmentColumn(), config.getSeqColumn(), config.getTimeColumn(), config.getTimePattern());
		if (tracks.isEmpty()) {
			String errorMsg = "No track found in the provided track shapefile.";
			LOGGER.severe(errorMsg);
			throw new RuntimeException(errorMsg);
		}
		LOGGER.info("Read " + tracks.size() + " tracks.");

		// add all tracks to manager
		TrackManager manager = new TrackManager();
		for (Track t : tracks)
			manager.addTrack(t);

		// ################## read the outer frame from shapefile ##################
		Polygon outerFrame;
		if (config.getFrameFile().isPresent()) {
			List<Feature> frames = FeatureReader
					.readFeatureFromShapefile(config.getFrameFile().get().getAbsoluteFile());
			if (frames.isEmpty()) {
				String errorMsg = "No frame found in the provided frame shapefile.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			} else if (frames.size() > 1 || frames.get(0).getGeometry().getNumGeometries() > 1) {
				LOGGER.info("Found more than one frame, only the first one is used.");
			} else if (!frames.get(0).getGeometryType().endsWith("Polygon")) {
				String errorMsg = "The provided frame is not a polygon!";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}

			outerFrame = (Polygon) frames.get(0).getGeometry().getGeometryN(0);
		} else {
			outerFrame = Tools.createEnvelope(config.getPseudomercatorFrameBounds().get());
		}

		// ################## generate the inner frame ##################
		double minExtent = outerFrame.getEnvelopeInternal().getHeight();
		if (outerFrame.getEnvelopeInternal().getWidth() < minExtent) {
			minExtent = outerFrame.getEnvelopeInternal().getWidth();
			LOGGER.fine("Minimum extent is width: " + minExtent + " [m]");
		} else {
			LOGGER.fine("Minimum extent is height: " + minExtent + " [m]");
		}

		// distance from outer to inner frame
		int b;
		if (config.getB().isPresent()) {
			b = config.getB().get();
		} else {
			b = (int) Math.ceil(config.getBeta().get() * minExtent);
		}
		LOGGER.fine("b = " + b + " [m]");

		Polygon innerFrame = (Polygon) outerFrame.buffer(-b);

		// ################## create times ##################
		double[][] lookbackUnits = new double[config.getNumShapes()][2];
		double lookbackMax = Double.NaN;
		LookBackType lookbackType = null;
		if (config.getMaxTime().isPresent()) {
			lookbackType = LookBackType.TIME;
			lookbackMax = config.getMaxTime().get();
			LOGGER.fine("lookbackMax = " + lookbackMax + " [s]");
		} else if (config.getMaxDist().isPresent()) {
			lookbackType = LookBackType.DISTANCE;
			lookbackMax = config.getMaxDist().get();
			LOGGER.fine("lookbackMax = " + lookbackMax + " [m]");
		} else {
			lookbackType = LookBackType.RANGE;
			lookbackMax = config.getMaxRange().get();
			LOGGER.fine("lookbackMax = " + lookbackMax + " [m] - range!");
		}
		lookbackUnits = generateLookbackUnits(lookbackMax, config.getNumShapes(), config.getOverlap());

		// ################## generate the scale factor ##################
		double scaleFactor;
		if (config.getScaleFactor().isPresent()) {
			scaleFactor = config.getScaleFactor().get();
		} else {
			scaleFactor = b / lookbackMax;
		}
		if (config.getMaxSpeed().isPresent()) {
			LOGGER.fine("max_speed = " + config.getMaxSpeed().get() + " [m/s] = " + config.getMaxSpeed().get() * 3.6
					+ " [km/h]");

			scaleFactor = b / (config.getMaxSpeed().get() * lookbackMax);
			if (config.getKappa().isPresent())
				scaleFactor /= config.getKappa().get();
		}

		LOGGER.fine("s = " + scaleFactor);

		// ################## all input arguments are parsed here ##################

		// ################## starting analysis for each track ##################
		ArrayList<Track> queriedTracks = manager.queryTracks(outerFrame.getEnvelopeInternal());
		LOGGER.info("Found " + queriedTracks.size() + " tracks inside outer frame.");

		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		LinkedList<Geometry> crossingpoints = new LinkedList<>();

		LinkedList<EvolutionSymbol> signatureList = new LinkedList<>();
		for (Track track : queriedTracks) {
			LOGGER.info("Starting track: " + track.getId());

			// get viewpoint crossings for track: intersections with inner and outer frame
			LOGGER.finest("Crossings: ");
			List<TrackPoint[]> pairs = track.getViewpointCrossings(outerFrame, innerFrame);
			for (TrackPoint[] pair : pairs)
				LOGGER.finest(pair[0] + " - " + pair[1]);
			LOGGER.fine("the track has " + pairs.size() + " intersections with the map window");
			// -> the first point of the pair is the inner point where the shape is computed

			// for each inner intersection point
			for (TrackPoint[] tp_pair : pairs) {

				// only consider complete pairs
				if (tp_pair[0] != null && tp_pair[1] != null) {

					TrackPoint current = tp_pair[0]; // here, cakeslice has to be computed
					LOGGER.finest("Current point: " + current);
					PointType tp_type = current.getType();

					TrackPoint previous = null;
					Direction direction = null;
					if (tp_type == PointType.CROSSING_FORWARD) {
						previous = track.getLatestPointBefore(current.getTimestamp());
						direction = Direction.FORWARD;
					} else if (tp_type == PointType.CROSSING_BACKWARD) {
						previous = track.getEarliestPointAfter(current.getTimestamp());
						direction = Direction.BACKWARD;
					} else {
						LOGGER.severe("Unexpected point type in pair.");
					}
					LOGGER.fine("cakeslice has to be computed " + direction);
					LOGGER.fine("Previous point: " + previous);

					LinkedList<EvolutionSymbol> signature = ComputeEvolution.computeEvolution(direction, track,
							lookbackUnits, config.getKappa().get(), current, previous, SignatureType.CAKE_SLICE,
							config.getVariant(), scaleFactor, lookbackType);
					signatureList.addAll(signature);

					Point crossingpoint = gf.createPoint(tp_pair[1].getCoordinate());
					crossingpoints.add(crossingpoint);
				}
			}
		}

		// ################## exort results ##################
		LinkedList<Double> endtimes = new LinkedList<>();
		for (int j = 0; j < (signatureList.size() / config.getNumShapes()); j++) {
			for (int i = 0; i < lookbackUnits.length; i++) {
				endtimes.add(lookbackUnits[i][1]);
			}
		}

		LinkedList<Integer> id = new LinkedList<>();
		for (int j = 0; j < (signatureList.size() / config.getNumShapes()); j++) {
			for (int i = 1; i <= config.getNumShapes(); i++) {
				id.add(i);
			}
		}

		EvolutionSymbol.writeSymbolWithPointsToShp(//
				new File(config.getOutputDir(), "offscreen"), //
				signatureList, //
				endtimes, //
				null);

		LinkedList<Geometry> frames = new LinkedList<>();
		frames.add(outerFrame);
		frames.add(innerFrame);
		frames.add(outerFrame.difference(innerFrame));
		LinkedList<Double> empty = new LinkedList<>();
		empty.add(0d);
		empty.add(0d);
		empty.add(0d);
		LinkedList<Integer> types = new LinkedList<>();
		types.add(0);
		types.add(1);
		types.add(2);
		Tools.exportGeometryAsPolygons(new File(config.getOutputDir(), "frames.shp").getAbsolutePath(), frames, empty,
				types);

		LOGGER.info("Offscreen-evaluation app finished successfully.");
	}

	public static double unitsPerShapeFromLookbackMax(double lookbackMax, int numShapes, double overlap) {
		return lookbackMax / (numShapes - overlap * (numShapes - 1));
	}

	public static double[][] generateLookbackUnits(double lookbackMax, int numShapes, double overlap) {
		double unitsPerShape = unitsPerShapeFromLookbackMax(lookbackMax, numShapes, overlap);
		double[][] lookbackUnits = new double[numShapes][2];
		lookbackUnits[0][0] = 0;
		lookbackUnits[0][1] = unitsPerShape;
		for (int i = 1; i < numShapes; i++) {
			lookbackUnits[i][0] = lookbackUnits[i - 1][0] + (1 - overlap) * unitsPerShape;
			lookbackUnits[i][1] = lookbackUnits[i][0] + unitsPerShape;
		}
		return lookbackUnits;
	}

	public static double lookbackMaxFromUnitsPerShape(double unitsPerShape, int numShapes, double overlap) {
		return unitsPerShape + (numShapes - 1) * (1 - overlap) * unitsPerShape;
	}

	public static double[][] generateLookbackUnitsFromUnitsPerShape(double unitsPerShape, int numShapes,
			double overlap) {
		double[][] lookbackUnits = new double[numShapes][2];
		lookbackUnits[0][0] = 0.0; // starts at current point
		lookbackUnits[0][1] = unitsPerShape; // end of first interval
		for (int i = 1; i < numShapes; i++) {
			lookbackUnits[i][0] = lookbackUnits[i - 1][1] - overlap * unitsPerShape;
			lookbackUnits[i][1] = lookbackUnits[i][0] + unitsPerShape;
		}
		return lookbackUnits;
	}
}
