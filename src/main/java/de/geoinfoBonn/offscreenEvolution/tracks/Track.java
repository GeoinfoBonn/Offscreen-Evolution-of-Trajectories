package de.geoinfoBonn.offscreenEvolution.tracks;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FeatureSource;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.linearref.LinearLocation;
import org.locationtech.jts.linearref.LocationIndexedLine;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import de.geoinfoBonn.offscreenEvolution.tracks.TrackPoint.PointType;

public final class Track {

	private static final Logger LOGGER = Logger.getLogger(Track.class.getName());

	private String id;
	private int segmentId;
	private LinkedList<TrackPoint> points;

	private LineString line;
	private LocationIndexedLine locationIndexedLine;
	private LengthIndexedLine lengthIndexedLine;

	private boolean finalized;

	/**
	 * Creates a new track.
	 * 
	 * @param id        unique id of the track
	 * @param segmentId counter for tracks that have several segments
	 */
	public Track(String id, int segmentId) {
		super();
		this.id = id;
		this.segmentId = segmentId;
		this.points = new LinkedList<>();
		this.finalized = false;
	}

	/**
	 * Creates a (non-finalized) copy of the given track.
	 * 
	 * @param track track to copy
	 */
	public Track(Track track) {
		this.id = track.id;
		this.segmentId = track.segmentId;
		this.points = new LinkedList<>(track.points);
		this.finalized = false;
	}

	/**
	 * Adds a new point to this track. Can only be done before track is finalized.
	 * 
	 * @param point new point to add
	 */
	public void addTrackPoint(TrackPoint point) {
		this.checkFinalized(false, "Point cannot be added when Track is already finalized.");

		this.points.add(point);
	}

	public List<TrackPoint> getPointsInRange(TrackPoint p, double range) {
		this.checkFinalized(true);

		boolean searchBefore = range < 0;
		range = Math.abs(range);

		LinkedList<TrackPoint> extractedPoints = new LinkedList<>();
		boolean foundStart = searchBefore;
		for (TrackPoint tp : points) {
			if (!foundStart)
				if (tp.getSequence() >= p.getSequence())
					foundStart = true;
				else
					continue;

			if (searchBefore && tp.getSequence() > p.getSequence())
				break;

			if (p.getCoordinate().distance(tp.getCoordinate()) <= range)
				extractedPoints.add(tp);
			else if (!searchBefore && foundStart)
				break;
			else if (searchBefore && !extractedPoints.isEmpty())
				extractedPoints = new LinkedList<>();

			if (searchBefore && tp.getSequence() == p.getSequence())
				break;
		}

		return extractedPoints;
	}

	/**
	 * Queries all points that are between the given point on the track in a
	 * distance of maximum {@code length}. Negative length looks backwards.
	 * 
	 * @param p      starting point
	 * @param length distance to look for (negative = backwards)
	 * @return all points in the specified interval
	 */
	public List<TrackPoint> getPointsBetween(TrackPoint p, double length) {
		this.checkFinalized(true);

		double lengthTp = lengthIndexedLine.indexOf(p.getCoordinate());
		double length2 = lengthTp + length;

		LineString extractedSection = null;
		if (lengthTp < length2)
			extractedSection = (LineString) lengthIndexedLine.extractLine(lengthTp, length2);
		else
			extractedSection = (LineString) lengthIndexedLine.extractLine(length2, lengthTp);

		LinearLocation start = locationIndexedLine.indexOf(extractedSection.getStartPoint().getCoordinate());
		LinearLocation end = locationIndexedLine.indexOf(extractedSection.getEndPoint().getCoordinate());

		double startSeq = start.getSegmentIndex() + start.getSegmentFraction();
		double endSeq = end.getSegmentIndex() + end.getSegmentFraction();

		LinkedList<TrackPoint> extractedPoints = new LinkedList<>();
		boolean foundStart = false;
		for (TrackPoint tp : points) {
			if (!foundStart)
				if (tp.getSequence() >= startSeq)
					foundStart = true;
				else
					continue;

			if (tp.getSequence() > endSeq)
				break;

			extractedPoints.add(tp);

			if (tp.getSequence() == endSeq)
				break;
		}

		return extractedPoints;
	}

	/**
	 * Queries all track points inside the specified time range.
	 * 
	 * @param startInclusive start of time range
	 * @param endExclusive   end of time range
	 * @return all track points in the specified time range
	 */
	public List<TrackPoint> getPointsBetween(LocalDateTime startInclusive, LocalDateTime endExclusive) {
		this.checkFinalized(true);

		LinkedList<TrackPoint> extractedPoints = new LinkedList<>();
		boolean foundStart = false;
		for (TrackPoint tp : points) {
			if (tp.getTimestamp().isAfter(endExclusive) || tp.getTimestamp().isEqual(endExclusive))
				break;

			if (!foundStart) {
				if (tp.getTimestamp().isAfter(startInclusive) || tp.getTimestamp().isEqual(startInclusive)) {
					foundStart = true;
					extractedPoints.add(tp);
				}
				continue;
			}

			extractedPoints.add(tp);
		}

		return extractedPoints;
	}

	/**
	 * Queries all points that are a certain time span before a given time.
	 * 
	 * @param endInclusive time
	 * @param timespan     time span
	 * @return all track points in the specified time range
	 */
	public List<TrackPoint> getPointsBefore(LocalDateTime endInclusive, Duration timespan) {
		return getPointsBetween(endInclusive.minus(timespan), endInclusive);
	}

	/**
	 * Queries all points that are a certain time span after a given time.
	 * 
	 * @param startInclusive time
	 * @param timespan       time span
	 * @return all track points in the specified time range
	 */
	public List<TrackPoint> getPointsAfter(LocalDateTime startInclusive, Duration timespan) {
		return getPointsBetween(startInclusive, startInclusive.plus(timespan));
	}

	public TrackPoint getLatestPointBefore(LocalDateTime time) {
		this.checkFinalized(true);

		TrackPoint previous = null;
		for (TrackPoint tp : points) {
			if (!tp.getTimestamp().isBefore(time))
				return previous;
			previous = tp;
		}
		return previous;
	}

	public TrackPoint getEarliestPointAfter(LocalDateTime time) {
		this.checkFinalized(true);

		for (TrackPoint tp : points) {
			if (tp.getTimestamp().isAfter(time))
				return tp;
		}
		return null;
	}

	/**
	 * Calculates all segments where this track crosses the inner and the outer
	 * envelope.
	 * 
	 * @param outer outer envelope
	 * @param inner inner envelope
	 * @return list of pairs with pair[0] --> inner crossing, pair[1] --> outer
	 *         crossing
	 */
	public List<TrackPoint[]> getViewpointCrossings(Envelope outer, Envelope inner) {
		Polygon outerP = Track.envelopeToPolygon(outer);
		Polygon innerP = Track.envelopeToPolygon(inner);

		return getViewpointCrossings(outerP, innerP);
	}

	/**
	 * Calculates all segments where this track crosses the inner and the outer
	 * envelope.
	 * 
	 * @param outer outer envelope
	 * @param inner inner envelope
	 * @return
	 */
	public List<TrackPoint[]> getViewpointCrossings(Polygon outer, Polygon inner) {
		this.checkFinalized(true);

		if (!outer.contains(inner))
			throw new IllegalArgumentException("Inner polygon must be strictly inside the outer polygon.");

		// collect all the crossing points
		List<TrackPoint> innerCrossings = this.calculateCrossingPoints(inner);
		List<TrackPoint> outerCrossings = this.calculateCrossingPoints(outer, true);

		List<TrackPoint> allCrossings = new LinkedList<>();
		allCrossings.addAll(innerCrossings);
		allCrossings.addAll(outerCrossings);
		Collections.sort(allCrossings);

		// initialize result
		List<TrackPoint[]> crossingPairs = new LinkedList<>();

		// special case: starting point is outside of inner but inside outer polygon
		boolean startIsInInner = inner.contains(line.getStartPoint());
		boolean startIsInOuter = outer.contains(line.getStartPoint());
		boolean startInBorderRegion = !startIsInInner && startIsInOuter;
		if (startInBorderRegion) {
//			System.err.println(
//					"Please check if true! Special case: starting point is outside of inner but inside outer polygon");
			TrackPoint[] pair = new TrackPoint[2];
			if (innerCrossings.contains(allCrossings.get(0)))
				pair[0] = allCrossings.get(0);
			else
				pair[1] = allCrossings.get(0);
			crossingPairs.add(pair);
			allCrossings.remove(0);
		}

		// special case: ending point is outside of inner but inside outer polygon
		if (allCrossings.size() % 2 == 1) {
//			System.err.println(
//					"Please check if true! Special case: ending point is outside of inner but inside outer polygon");
			TrackPoint[] pair = new TrackPoint[2];
			if (innerCrossings.contains(allCrossings.get(allCrossings.size() - 1)))
				pair[0] = allCrossings.get(allCrossings.size() - 1);
			else
				pair[1] = allCrossings.get(allCrossings.size() - 1);
			crossingPairs.add(pair);
			allCrossings.remove(allCrossings.size() - 1);
		}

		List<TrackPoint> toRemove = new LinkedList<>();
		PointType expectedType = null;
		for (int i = 0; i < allCrossings.size(); ++i) {
			if (expectedType == null) {
//				expectedType = allCrossings.get(i).getType();
				expectedType = allCrossings.get(i).getType() == PointType.CROSSING_OUTER ? PointType.CROSSING_BACKWARD
						: PointType.CROSSING_OUTER;
			} else {
				if (expectedType != allCrossings.get(i).getType()) {
					toRemove.add(allCrossings.get(i - 1));
					toRemove.add(allCrossings.get(i));
				}
				expectedType = null;
			}
		}
		allCrossings.removeAll(toRemove);

		// transform crossing points to pairs of crossings
		for (int i = 0; i < allCrossings.size(); i = i + 2) {
			TrackPoint[] pair = new TrackPoint[2];
			if (allCrossings.get(i).getType() == PointType.CROSSING_FORWARD) {
				pair[0] = allCrossings.get(i);
				pair[1] = allCrossings.get(i + 1);
			} else {
				pair[0] = allCrossings.get(i + 1);
				pair[1] = allCrossings.get(i);
			}
			crossingPairs.add(pair);
		}
		return crossingPairs;
	}

	public static Polygon envelopeToPolygon(Envelope env) {
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		Coordinate topLeft = new Coordinate(env.getMinX(), env.getMaxY());
		Coordinate topRight = new Coordinate(env.getMaxX(), env.getMaxY());
		Coordinate bottomLeft = new Coordinate(env.getMinX(), env.getMinY());
		Coordinate bottomRight = new Coordinate(env.getMaxX(), env.getMinY());
		return gf.createPolygon(new Coordinate[] { topLeft, topRight, bottomRight, bottomLeft, topLeft });
	}

	public List<TrackPoint> calculateCrossingPoints(Envelope env) {
		return calculateCrossingPoints(env, false);
	}

	public List<TrackPoint> calculateCrossingPoints(Envelope env, boolean isOuter) {
		this.checkFinalized(true);

		return calculateCrossingPoints(Track.envelopeToPolygon(env), isOuter);
	}

	public List<TrackPoint> calculateCrossingPoints(Polygon polygon) {
		return calculateCrossingPoints(polygon, false);
	}

	/**
	 * Calculates all crossing points of this Track with the given polygon.
	 * 
	 * @param polygon polygon to calculate crossing points with
	 * @return all found crossings
	 */
	public List<TrackPoint> calculateCrossingPoints(Polygon polygon, boolean isOuter) {
		this.checkFinalized(true);

		LinkedList<TrackPoint> crossings = new LinkedList<>();

		Geometry intersections = line.intersection(polygon.getBoundary());
		if (intersections.isEmpty())
			return crossings;
		for (int i = 0; i < intersections.getNumGeometries(); ++i) {
			if (!intersections.getGeometryN(i).getGeometryType().equals("Point")) {
				LOGGER.warning("Why is this the case?");
				continue;
			}
			Point intersection = (Point) intersections.getGeometryN(i);

			LinearLocation l = locationIndexedLine.indexOf(intersection.getCoordinate());

			double fraction = l.getSegmentFraction();
			Coordinate location = l.getCoordinate(line);
			TrackPoint prev = points.get(l.getSegmentIndex());
			TrackPoint next = points.size() > l.getSegmentIndex() + 1 ? points.get(l.getSegmentIndex() + 1) : null;

			LocalDateTime interpolatedTime = null;
			if (prev.getTimestamp() != null && next != null && next.getTimestamp() != null) {
				long nanosBetween = Duration.between(prev.getTimestamp(), next.getTimestamp()).toNanos();
				interpolatedTime = prev.getTimestamp().plus(Duration.ofNanos((long) (nanosBetween * fraction)));
			} else if (next == null) {
				interpolatedTime = prev.getTimestamp();
			}

			TrackPoint tp = new TrackPoint(location, interpolatedTime,
					points.get(l.getSegmentIndex()).getSequence() + fraction);
			crossings.add(tp);
		}

		Collections.sort(crossings);
		boolean nextIsForward = polygon.contains(line.getStartPoint());
		for (TrackPoint p : crossings) {
			PointType type = nextIsForward ? PointType.CROSSING_FORWARD : PointType.CROSSING_BACKWARD;
			if (isOuter)
				type = PointType.CROSSING_OUTER;
			else
				nextIsForward = !nextIsForward; // initialize next iteration

			p.setType(type);
		}

		return crossings;
	}

	/**
	 * Finalizes this line and creates all relevant indices for efficient queries.
	 */
	public void finalizeLine() {
		Collections.sort(points);
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();

		Coordinate[] c = new Coordinate[points.size()];
		for (int i = 0; i < points.size(); ++i) {
			c[i] = points.get(i).getCoordinate();
		}
//		this.minSequence = points.getFirst().getSequence();
		this.line = gf.createLineString(c);
		this.locationIndexedLine = new LocationIndexedLine(line);
		this.lengthIndexedLine = new LengthIndexedLine(line);
		this.finalized = true;
	}

	public static List<Track> importFromShapefile(String filename) {
		return importFromShapefile(filename, null, null, null, null, null);
	}

	public String getId() {
		return id;
	}

	public int getSegmentId() {
		return segmentId;
	}

	public LinkedList<TrackPoint> getPoints() {
		return points;
	}

	@Override
	public String toString() {
		return "Track [id=" + id + ", segmentId=" + segmentId + ", points=" + points.size() + "]";
	}

	public LineString getAsLineString() {
		return line;
	}

	public static List<Track> importFromShapefile(String filename, String idColumn, String subtrackIdColumn,
			String seqColumn, String timeColumn, String timePattern) {
		if (!filename.endsWith(".shp")) {
			System.out.println("Provides file is not a shapefile (filename does not end with .shp)");
			return null;
		}

		HashMap<String, Track> tracks = new HashMap<>();
		try {

			File shpfile = new File(filename);
			Map<String, Object> map = new HashMap<>();
			map.put("url", shpfile.toURI().toURL());

			DataStore dataStore = DataStoreFinder.getDataStore(map);
			String typeName = dataStore.getTypeNames()[0];

			FeatureSource<SimpleFeatureType, SimpleFeature> source = dataStore.getFeatureSource(typeName);
			Filter filter = Filter.INCLUDE; // ECQL.toFilter("BBOX(THE_GEOM, 10,20,30,40)")

			FeatureCollection<SimpleFeatureType, SimpleFeature> myFeatureCollection = source.getFeatures(filter);
			try (FeatureIterator<SimpleFeature> features = myFeatureCollection.features()) {
				while (features.hasNext()) {
					SimpleFeature myFeature = features.next();
					String trackId = idColumn != null ? "" + myFeature.getAttribute(idColumn) : myFeature.getID();
					Geometry myGeometry = (Geometry) myFeature.getDefaultGeometry();

					if (myGeometry.getGeometryType().endsWith("LineString")) {
						for (int geoIndex = 0; geoIndex < myGeometry.getNumGeometries(); geoIndex++) {
							if (tracks.containsKey(id(trackId, geoIndex)))
								throw new IllegalArgumentException("ID already present");
							Track track = new Track(trackId, geoIndex);

							LineString geometry = (LineString) myGeometry.getGeometryN(geoIndex);
							Coordinate[] xyz = geometry.getCoordinates();
							for (int j = 0; j < xyz.length; ++j) {
								track.addTrackPoint(new TrackPoint(xyz[j], null, j));
							}

							tracks.put(id(trackId, geoIndex), track);
						}
					} else if (myGeometry.getGeometryType().equals("Point")) {
						int subtrackId = (int) (long) myFeature.getAttribute(subtrackIdColumn);
						if (!tracks.containsKey(id(trackId, subtrackId)))
							tracks.put(id(trackId, subtrackId), new Track(trackId, subtrackId));
						Track track = tracks.get(id(trackId, subtrackId));

						Point point = (Point) myGeometry;
						Coordinate xyz = point.getCoordinate();

						int point_seq = (int) (long) myFeature.getAttribute(seqColumn);
						LocalDateTime timestamp = LocalDateTime.parse((String) myFeature.getAttribute(timeColumn),
								DateTimeFormatter.ofPattern(timePattern));

						track.addTrackPoint(new TrackPoint(xyz, timestamp, point_seq));

					} else {
						System.err.println("The shapefile contains an unimplemented geometry type: "
								+ myGeometry.getGeometryType());
					}
				}
			} catch (Exception ex) {
				throw new RuntimeException(ex);
			} finally {
				dataStore.dispose();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		List<Track> result = tracks.values().stream().collect(Collectors.toList());
		for (Track t : result)
			t.finalizeLine();
		return result;
	}

	public static void writeToShapeFile(String filename, List<Track> tracks, String idColumn, String subtrackIdColumn,
			String seqColumn, String timeColumn, String timePattern, CoordinateReferenceSystem crs) {

		long currentTime = System.currentTimeMillis();
		SimpleFeatureType featureType = createFeatureType(idColumn, subtrackIdColumn, seqColumn, timeColumn, crs);

		DefaultFeatureCollection collection = new DefaultFeatureCollection();
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory(null);
		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		for (Track track : tracks) {
			for (TrackPoint tp : track.points) {
				Point p = gf.createPoint(tp.getCoordinate());

				featureBuilder.add(p);
				featureBuilder.add(track.getId());
				featureBuilder.add(track.getSegmentId());
				featureBuilder.add(tp.getSequence());
				featureBuilder.add(tp.getTimestamp().format(DateTimeFormatter.ofPattern(timePattern)));

				SimpleFeature feature = featureBuilder.buildFeature(null);
				collection.add(feature);
			}
		}

		try {
			File output = new File(filename);
			if (output.getParentFile() != null) {
				output.getParentFile().mkdirs();
			}
			output.createNewFile();
			ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();

			Map<String, Serializable> params = new HashMap<String, Serializable>();
			params.put("url", output.toURI().toURL());
			params.put("create spatial index", Boolean.TRUE);

			ShapefileDataStore newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);
			newDataStore.createSchema(featureType);

			/*
			 * Write the features to the shapefile
			 */
			Transaction transaction = new DefaultTransaction("create");

			String typeName = newDataStore.getTypeNames()[0];
			SimpleFeatureSource featureSource = newDataStore.getFeatureSource(typeName);

			if (featureSource instanceof SimpleFeatureStore) {
				SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;

				featureStore.setTransaction(transaction);
				try {
					featureStore.addFeatures(collection);
					transaction.commit();

				} catch (Exception problem) {
					problem.printStackTrace();
					transaction.rollback();

				} finally {
					transaction.close();
				}
			} else {
				System.out.println(typeName + " does not support read/write access");
			}

			System.out.println("Shape written to " + filename + " in "
					+ (System.currentTimeMillis() - currentTime) / 1000.0 + "s.");
		} catch (MalformedURLException e) {
			System.err.println("Error writing shape.");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error writing shape.");
			e.printStackTrace();
		}
	}

	private static SimpleFeatureType createFeatureType(String idColumn, String subtrackIdColumn, String seqColumn,
			String timeColumn, CoordinateReferenceSystem crs) {

		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName("Location");
		if (crs == null)
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		else
			builder.setCRS(crs); // <- Coordinate reference system

		// add attributes in order
		builder.add("the_geom", Point.class);
		builder.add(idColumn, String.class);
		builder.add(subtrackIdColumn, Integer.class);
		builder.add(seqColumn, Integer.class);
		builder.add(timeColumn, String.class);

		// build the type
		final SimpleFeatureType LOCATION = builder.buildFeatureType();

		return LOCATION;
	}

	public static String id(String trackId, int subtrackId) {
		return trackId + "_" + String.format("%03d", subtrackId);
	}

	public void checkFinalized(boolean expectedState) {
		String message = expectedState ? "Finalize before using this function!"
				: "Function can only be used before finalizing";
		checkFinalized(expectedState, message);
	}

	public void checkFinalized(boolean expectedState, String message) {
		if (finalized != expectedState)
			throw new IllegalStateException(message);
	}
}
