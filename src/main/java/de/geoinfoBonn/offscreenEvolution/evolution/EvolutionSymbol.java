package de.geoinfoBonn.offscreenEvolution.evolution;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.geotools.data.shapefile.ShapefileDumper;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import de.geoinfoBonn.offscreenEvolution.utils.GeometricTransformations;
import de.geoinfoBonn.offscreenEvolution.utils.Tools;

public abstract class EvolutionSymbol {
	private static final Logger LOGGER = Logger.getLogger(EvolutionSymbol.class.getName());

	static final AtomicLong NEXT_ID = new AtomicLong(0);
	final long id = NEXT_ID.getAndIncrement(); // unique id of the symbol
	Geometry symbol; // shape of the symbol
	ArrayList<Point> relevantPoints;
	LinkedList<Integer> coveredPoints;

	public static class CakesliceSymbol extends EvolutionSymbol {

		private static final double EPS_RANGE = 1e-9;
		private static final double EPS_ANGLE = 1e-12;

		private Coordinate sourceLocation;
		private double openeningAngle;
		private double range;
		private double directionOfRightEdge;

		public CakesliceSymbol(Coordinate sourceLocation, double openingAngle, double range,
				double directionOfRightEdge, double scaleFactor, List<Point> relevantPoints) {
			// check input
			if (openingAngle <= 0 || range <= 0)
				throw new IllegalArgumentException("Opening angle and range must be larger than zero.");
			if (directionOfRightEdge < 0 || directionOfRightEdge >= 2 * Math.PI)
				throw new IllegalArgumentException("Direction of right edge must be in [0, 2*pi)");
			if (scaleFactor <= 0 || scaleFactor > 1)
				throw new IllegalArgumentException("Scale factor must be in (0, 1]");

			// set given values
			this.sourceLocation = sourceLocation;
			this.openeningAngle = openingAngle;
			this.range = range;
			this.directionOfRightEdge = directionOfRightEdge;

			// compute derived values
			this.symbol = Tools.createCakesliceSymbol(directionOfRightEdge, sourceLocation, range, openingAngle,
					scaleFactor);
			this.setRelevantPoints(relevantPoints);
			this.computeCoveredPoints();

			LOGGER.finest(n() + "," + k());
		}

		private void computeCoveredPoints() {
			this.coveredPoints = new LinkedList<>();
			for (int i = 0; i < n(); ++i) {
				Coordinate c = relevantPoints.get(i).getCoordinate();
				if (this.covers(c))
					coveredPoints.add(i);
			}
		}

		public boolean covers(Coordinate c) {
			double r = sourceLocation.distance(c);
			if (r > range + EPS_RANGE)
				return false;

			double direction = GeometricTransformations.angleFromX(sourceLocation, c);
			if (direction < directionOfRightEdge - EPS_ANGLE
					|| direction > directionOfRightEdge + openeningAngle + EPS_ANGLE)
				return false;

			return true;
		}

	}

	public void setRelevantPoints(List<Point> relevantPoints) {
		Objects.requireNonNull(relevantPoints);
		this.relevantPoints = new ArrayList<>(relevantPoints);
	}

	public long getId() {
		return id;
	}

	public double getPerimeter() {
		return symbol.getBoundary().getLength();
	}

	public double getArea() {
		return symbol.getArea();
	}

	public Geometry getSymbol() {
		return symbol;
	}

	public List<Point> getRelevantPoints() {
		return relevantPoints;
	}

	public List<Point> getCoveredPoints() {
		return coveredPoints.stream().map(i -> relevantPoints.get(i)).collect(Collectors.toUnmodifiableList());
	}

	public LinkedList<Integer> getCoveredPointIds() {
		return coveredPoints;
	}

	public int n() {
		return relevantPoints.size();
	}

	public int k() {
		return coveredPoints.size();
	}

	public static void writeSymbolWithPointsToShp(File outputFile, List<EvolutionSymbol> symbols,
			LinkedList<Double> time, CoordinateReferenceSystem crs) {
		if (crs == null)
			try {
				crs = CRS.decode("EPSG:25832");
			} catch (NoSuchAuthorityCodeException e) {
				e.printStackTrace();
			} catch (FactoryException e) {
				e.printStackTrace();
			}

		long starttime = System.currentTimeMillis();
		SimpleFeatureType sFeatureType = EvolutionSymbol.createSymbolFeatureType(crs, time != null);
		SimpleFeatureType pFeatureType = EvolutionSymbol.createPointsFeatureType(crs, time != null);
		DefaultFeatureCollection sCollection = EvolutionSymbol.createSymbolFeatureCollection(sFeatureType, symbols,
				time);
		DefaultFeatureCollection pCollection = EvolutionSymbol.createPointsFeatureCollection(pFeatureType, symbols,
				time);
		EvolutionSymbol.writeToFile(outputFile.getParentFile(), outputFile.getName() + "_symbols", sCollection);
		EvolutionSymbol.writeToFile(outputFile.getParentFile(), outputFile.getName() + "_points", pCollection);

		LOGGER.info("Symbols and points written to " + outputFile + " in "
				+ (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
	}

	public static void writeSymbolsToShp(File outputFile, ArrayList<EvolutionSymbol> symbols, LinkedList<Double> time,
			CoordinateReferenceSystem crs) {
		long starttime = System.currentTimeMillis();
		SimpleFeatureType sFeatureType = EvolutionSymbol.createSymbolFeatureType(crs, time != null);
		DefaultFeatureCollection sCollection = EvolutionSymbol.createSymbolFeatureCollection(sFeatureType, symbols,
				time);
		EvolutionSymbol.writeToFile(outputFile.getParentFile(), outputFile.getName(), sCollection);

		LOGGER.info(
				"Symbols written to " + outputFile + " in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
	}

	private static SimpleFeatureType createSymbolFeatureType(CoordinateReferenceSystem crs, boolean withTime) {
		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName("EvolutionSymbolBuilder");
		if (crs == null)
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		else
			builder.setCRS(crs); // <- Coordinate reference system

		// add attributes in order
		builder.add("the_geom", Polygon.class);
		builder.add("symbol_id", Integer.class);
		if (withTime)
			builder.add("time", Double.class);

		// build the type
		final SimpleFeatureType SFT = builder.buildFeatureType();
		return SFT;
	}

	private static SimpleFeatureType createPointsFeatureType(CoordinateReferenceSystem crs, boolean withTime) {
		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName("EvolutionPointsBuilder");
		if (crs == null)
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		else
			builder.setCRS(crs); // <- Coordinate reference system

		// add attributes in order
		builder.add("the_geom", Point.class);
		builder.add("symbol_id", Integer.class);
		builder.add("covered", Boolean.class);
		if (withTime)
			builder.add("time", Double.class);

		// build the type
		final SimpleFeatureType SFT = builder.buildFeatureType();
		return SFT;
	}

	private static DefaultFeatureCollection createSymbolFeatureCollection(SimpleFeatureType featureType,
			List<EvolutionSymbol> symbols, LinkedList<Double> time) {

		DefaultFeatureCollection collection = new DefaultFeatureCollection();
		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		int i = 0;
		for (EvolutionSymbol symbol : symbols) {
			SimpleFeature feature;

			featureBuilder.add(symbol.getSymbol());
			featureBuilder.add(symbol.getId());
			if (time != null)
				featureBuilder.add(time.get(i));
			feature = featureBuilder.buildFeature(null);
			collection.add(feature); // add symbol
			i++;
		}

		return collection;
	}

	private static DefaultFeatureCollection createPointsFeatureCollection(SimpleFeatureType featureType,
			List<EvolutionSymbol> symbols, LinkedList<Double> time) {

		DefaultFeatureCollection collection = new DefaultFeatureCollection();
		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		int j = 0;
		for (EvolutionSymbol symbol : symbols) {
			List<Point> relevantPoints = symbol.getRelevantPoints();
			List<Integer> coveredPoints = symbol.getCoveredPointIds();
			for (int i = 0; i < symbol.n(); ++i) {
				featureBuilder.add(relevantPoints.get(i));
				featureBuilder.add(symbol.getId());
				featureBuilder.add(coveredPoints.contains(i));
				if (time != null)
					featureBuilder.add(time.get(j));
				SimpleFeature feature = featureBuilder.buildFeature(null);
				collection.add(feature); // add relevant points
			}
			j++;
		}

		return collection;
	}

	private static void writeToFile(File outputDir, String filename, DefaultFeatureCollection collection) {
		try {
			ShapefileDumper dumper = new ShapefileDumper(outputDir);
			// optiona, set a target charset (ISO-8859-1 is the default)
			dumper.setCharset(Charset.forName("ISO-8859-15"));
			// split when shp or dbf reaches 100MB
			int maxSize = 100 * 1024 * 1024;
			dumper.setMaxDbfSize(maxSize);
			// actually dump data
			dumper.dump(filename, collection);
		} catch (IOException e) {
			LOGGER.severe("Error writing shape.");
			e.printStackTrace();
		}
	}
}
