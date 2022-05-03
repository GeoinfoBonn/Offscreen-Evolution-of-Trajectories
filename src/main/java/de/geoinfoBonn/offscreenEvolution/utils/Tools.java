package de.geoinfoBonn.offscreenEvolution.utils;

import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Logger;

import org.geotools.data.DefaultTransaction;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geojson.feature.FeatureJSON;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.util.GeometricShapeFactory;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import de.geoinfoBonn.offscreenEvolution.io.Feature;

public class Tools {
	private static final Logger LOGGER = Logger.getLogger(Tools.class.getName());

	public static Polygon createCakesliceSymbol(double gamma, Point2D currentPoint, double r, double phi, double s) {
		return createCakesliceSymbol(gamma, new Coordinate(currentPoint.getX(), currentPoint.getY()), r, phi, s);
	}

	public static Polygon createCakesliceSymbol(double gamma, Coordinate currentPoint, double r, double phi, double s) {
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		GeometricShapeFactory shapeFactory = new GeometricShapeFactory(gf);
		shapeFactory.setCentre(currentPoint);

		LOGGER.finest("scaling factor for signature = " + s);
		shapeFactory.setSize(2 * r * s);
		LOGGER.finest("phi für cakeslice: " + phi * 180 / Math.PI);
		LineString arc = null;
		// if no opening angle due to n=k=1
		if (phi == 0.0) {
			System.err.println("Angle of cakeslice is 0.0 - please check");
			arc = shapeFactory.createArc(gamma, 0.001);
		}
		if (phi < 0) {
			arc = shapeFactory.createArc(gamma + phi, -phi);

		} else {
			arc = shapeFactory.createArc(gamma, phi); // start angle and size of angle
		}
		Coordinate[] coords = new Coordinate[arc.getCoordinates().length + 2];
		coords[0] = new Coordinate(currentPoint.getX(), currentPoint.getY());
		for (int i = 0; i < arc.getCoordinates().length; ++i) {
			coords[i + 1] = arc.getCoordinateN(i);
		}
		coords[coords.length - 1] = new Coordinate(currentPoint.getX(), currentPoint.getY());

		return gf.createPolygon(coords);
	}

	public static void exportAsGeoJSON(File file, List<Geometry> geometries, List<Double> time,
			CoordinateReferenceSystem inputCrs, CoordinateReferenceSystem outputCrs) {
		MathTransform transform = null;
		try {
			if (inputCrs == null)
				inputCrs = CRS.decode("EPSG:25832");
			if (outputCrs == null)
				outputCrs = CRS.decode("EPSG:3857");

			transform = CRS.findMathTransform(inputCrs, outputCrs);

		} catch (NoSuchAuthorityCodeException e1) {
			e1.printStackTrace();
		} catch (FactoryException e1) {
			e1.printStackTrace();
		}

		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(createFeatureType(outputCrs));
		SimpleFeatureBuilder pointFeatureBuilder = new SimpleFeatureBuilder(createCrossingFeatureType(outputCrs));

		FeatureJSON fjson = new FeatureJSON();
		DefaultFeatureCollection collection = new DefaultFeatureCollection();

//		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		SimpleFeature feature = null;
		int i = 0;
		for (Geometry geom : geometries) {
			Geometry geomTrans = null;
			try {
				geomTrans = JTS.transform(geom, transform);
			} catch (MismatchedDimensionException | TransformException e) {
				e.printStackTrace();
			}
			if (geom.getGeometryType().endsWith("Polygon")) {
				featureBuilder.add(geomTrans);
				featureBuilder.add(time.get(i));
				featureBuilder.add(0); // id
				feature = featureBuilder.buildFeature(null);
			} else if (geom instanceof Point) {
				pointFeatureBuilder.add(geomTrans);
				feature = pointFeatureBuilder.buildFeature(null);
			}

			collection.add(feature);
			i++;
		}

		try (BufferedWriter fileWriter = new BufferedWriter(new FileWriter(file))) {
			fjson.writeFeatureCollection(collection, fileWriter);
//			fileWriter.append(writer.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void exportGeometryAsPolygons(String filename, LinkedList<Geometry> geometries) {
		exportGeometryAsPolygons(filename, geometries, null, null, null);
	}

	public static void exportGeometryAsPolygons(String filename, LinkedList<Geometry> geometries,
			LinkedList<Double> time, LinkedList<Integer> id) {
		exportGeometryAsPolygons(filename, geometries, time, id, null);
	}

	public static void exportGeometryAsPolygons(String filename, LinkedList<Geometry> geometries,
			LinkedList<Double> time, LinkedList<Integer> id, CoordinateReferenceSystem crs) {

		Objects.requireNonNull(filename);
		Objects.requireNonNull(geometries);

		if (time != null) {
			if (geometries.size() != time.size()) {
				throw new IllegalArgumentException("Lengths of both lists have to be identical.");
			}
		}

		if (id != null) {
			if (geometries.size() != id.size()) {
				throw new IllegalArgumentException("Lengths of all three lists have to be identical.");
			}
		}

		long currentTime = System.currentTimeMillis();
		SimpleFeatureType featureType = null;
		try {
			if (crs == null)
				crs = CRS.decode("EPSG:25832");
			featureType = createFeatureType(crs);
		} catch (NoSuchAuthorityCodeException e1) {
			e1.printStackTrace();
		} catch (FactoryException e1) {
			e1.printStackTrace();
		}

		DefaultFeatureCollection collection = new DefaultFeatureCollection();

		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		int i = 0;
		for (Geometry geom : geometries) {
			if (geom instanceof LineString)
				featureBuilder.add(gf.createPolygon(((LineString) geom).getCoordinateSequence()));
			else if (geom.getGeometryType().endsWith("Polygon"))
				featureBuilder.add(geom);
			if (time != null) {
				featureBuilder.add(time.get(i));
			}
			if (id != null) {
				featureBuilder.add(id.get(i++));
			} else {
				featureBuilder.add(0);
				featureBuilder.add(0);
			}

			SimpleFeature feature = featureBuilder.buildFeature(null);
			collection.add(feature);
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
					newDataStore.dispose();
				}
			} else {
				LOGGER.severe(typeName + " does not support read/write access");
			}

			LOGGER.info("Shape written to " + filename + " in " + (System.currentTimeMillis() - currentTime) / 1000.0
					+ "s.");
		} catch (MalformedURLException e) {
			LOGGER.severe("Error writing shape.");
			e.printStackTrace();
		} catch (IOException e) {
			LOGGER.severe("Error writing shape.");
			e.printStackTrace();
		}
	}

	public static void exportLineStrings(String filename, LinkedList<LineString> ls) {
		long currentTime = System.currentTimeMillis();
		SimpleFeatureType featureType = null;
		try {
			featureType = createFeatureType(CRS.decode("EPSG:25832"));
		} catch (NoSuchAuthorityCodeException e1) {
			e1.printStackTrace();
		} catch (FactoryException e1) {
			e1.printStackTrace();
		}

		DefaultFeatureCollection collection = new DefaultFeatureCollection();

		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();

		for (LineString l : ls) {
			featureBuilder.add(gf.createPolygon(l.getCoordinateSequence()));
			featureBuilder.add(0);

			SimpleFeature feature = featureBuilder.buildFeature(null);
			collection.add(feature);
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
					newDataStore.dispose();
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

	private static SimpleFeatureType createFeatureType(CoordinateReferenceSystem crs) {

		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName("Location");
		if (crs == null)
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		else
			builder.setCRS(crs); // <- Coordinate reference system

		// add attributes in order
		builder.add("the_geom", MultiPolygon.class);
		builder.add("time", Double.class);
		builder.add("id", Integer.class);

		// build the type
		final SimpleFeatureType LOCATION = builder.buildFeatureType();

		return LOCATION;
	}

	private static SimpleFeatureType createCrossingFeatureType(CoordinateReferenceSystem crs) {

		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName("crosspoint");
		if (crs == null)
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		else
			builder.setCRS(crs); // <- Coordinate reference system

		// add attributes in order
		builder.add("the_geom", Point.class);

		// build the type
		final SimpleFeatureType CROSSPOINT = builder.buildFeatureType();

		return CROSSPOINT;
	}

	public static double getAverage(double[] marks) {
		double total = 0;
		for (int i = 0; i < marks.length; i++) {
			total += marks[i];
		}
		return total / marks.length;
	}

	// for median computation
	public static void swap(double[] a, int i1, int i2) {
		double temp = a[i1];
		a[i1] = a[i2];
		a[i2] = temp;
	}

	public static double median(double[] a, int from, int to) {
		int low = from;
		int high = to;
		int median = (low + high) / 2;
		do {
			if (high <= low) {
				return a[median];
			}
			if (high == low + 1) {
				if (a[low] > a[high]) {
					swap(a, low, high);
				}
				return a[median];
			}
			int middle = (low + high) / 2;
			if (a[middle] > a[high]) {
				swap(a, middle, high);
			}
			if (a[low] > a[high]) {
				swap(a, low, high);
			}
			if (a[middle] > a[low]) {
				swap(a, middle, low);
			}
			swap(a, middle, low + 1);
			int ll = low + 1;
			int hh = high;
			do {
				do {
					ll++;
				} while (a[low] > a[ll]);
				do {
					hh--;
				} while (a[hh] > a[low]);
				if (hh < ll) {
					break;
				}
				swap(a, ll, hh);
			} while (true);
			swap(a, low, hh);
			if (hh <= median) {
				low = ll;
			}
			if (hh >= median) {
				high = hh - 1;
			}
		} while (true);
	}

	// function must be adjusted to cakeslice shape!
	/**
	 * 
	 * @param k
	 * @param currentPoint
	 * @param sites
	 * @param gamma
	 * @param phi          result_tri[1]
	 * @param r            result_tri[0]
	 * @return local polar points that lie in cakeslice
	 */
	public static Point2D[] getPointsInCakeslice(int k, Point2D currentPoint, Point2D[] sites, double gamma, double phi,
			double r) {
		// transform currentPoint to local polar coordinates
		currentPoint = GeometricTransformations.global2local(gamma, currentPoint.getX(), currentPoint.getY(),
				currentPoint);
		currentPoint = GeometricTransformations.cart2pol(currentPoint); // FALSCH
		// store k points
		Point2D[] kpts = new Point2D[k];
		// iterate over all n points (sites)
		for (int i = 0; i < sites.length; i++) {
			Point2D p = sites[i];
			double alpha = Math.atan2(p.getY() - currentPoint.getY(), p.getX() - currentPoint.getX());
			if (alpha < 0) {
				alpha += 2 * Math.PI;
			}
			double dAlpha = Math.abs(alpha - gamma);
			double dist = currentPoint.distance(p);
			// check orientation
//			System.out.println("dAlpha = " + dAlpha + ", maxDif = " + maxDif);
			// only compare up to the 9. decimal place due to numerical accuracy
			double treshold = 0.0000000001;
			double currDif = (dAlpha - phi);
			double currDif2 = ((2 * Math.PI - dAlpha) - phi);
//			System.out.println("currDif = " + currDif + ", currDif2 = " + currDif2);
			if ((currDif <= treshold || currDif2 <= treshold)) {
				// if ((dAlpha <= maxDif || 2 * Math.PI - dAlpha <= maxDif)) { // point ON
				// triangle should also lie in!
				// check extent -> due to cakeslice: simple comparison with r
				double dist_gerechnet = r;
				double distDif = (dist - dist_gerechnet);
//				System.out.println("dist = " + dist);
//				System.out.println("dist_gerechnet = " + dist_gerechnet);
//				System.out.println("distDif: " + distDif);
				if ((distDif <= treshold) && dist > treshold) {
					// if (dist <= dist_gerechnet && dist != 0) { // current point must not lie in
//					counter++;
					kpts[i] = sites[i];
					System.out.println("POINT LIES IN");
				} else {
					System.err.println("Point does not count.");
				}
			}
		}
		return kpts;
	}

	// KLAPPT SO NICHT! EHER VERGLEICH MIT Lage bezï¿½glich KANTEN DES CAKESLICE?
	public static Point2D[] getPointsInCakesliceUsingCoords(int k, Point2D[] points, Coordinate[] cakeslicecoords) {
		Point2D[] kpts = new Point2D[k + 100]; // es sind genau k Punkte enthalten (testweise plus 100)
		int count = 0;
		for (int i = 0; i < points.length; i++) { // vergleiche alle points...
			Point2D p = points[i];
			for (int j = 0; j < cakeslicecoords.length; j++) { // ... mit allen cakeslicecoords
				Point2D c = new Point2D.Double(cakeslicecoords[j].x, cakeslicecoords[j].y);
				if (((p.getX() <= c.getX()) && (p.getY() <= c.getY()))
						|| ((p.getX() <= c.getX()) && (p.getY() >= c.getY()))
						|| ((p.getX() >= c.getX()) && (p.getY() <= c.getY()))
						|| ((p.getX() >= c.getX()) && (p.getY() >= c.getY()))) { // wenn point innerhalb der coords
																					// liegt, fï¿½ge in ktps hinzu
					kpts[count] = p;
					break;
				}
			}
			count++;
		}
		return kpts;
	}

	/**
	 * 
	 * @param ls
	 * @return
	 */
	public static List<Point2D> lineString2List(LineString ls) {
		List<Point2D> pts = new ArrayList<Point2D>();
		for (int i = 0; i < ls.getNumPoints(); i++) {
			pts.add(new Point2D.Double(ls.getCoordinateN(i).x, ls.getCoordinateN(i).y));

		}
		return pts;
	}

	/**
	 * 
	 * @param f imported shp feature: the GPS-track
	 * @return one linestring
	 */
	public static LineString preprocessLineString(Feature f) {
		// tracks to be processed
		LineString ls = null;

		// test amount of geometries
		int g = f.getGeometry().getNumGeometries();
		if (g > 1) {
			System.err.println("Only the first geometry of " + g + " geometries is processed");
		}

		// get first LineString (i.e. in case of MultiLineString)
		if (f.getGeometryType().endsWith("LineString")) {
			ls = (LineString) f.getGeometry().getGeometryN(0);
			System.out.println("# Points of track: " + ls.getNumPoints());
			System.out.println("Length of tracK [m]: " + ls.getLength());
		} else {
			System.err.println("Feature contains no LineString - continue");
			// continue;
		}
		return ls;
	}

	/**
	 * This function returns an index: at this position the new point has to be
	 * inserted into the linestring. It finds the correct position by comparing the
	 * distances between the existing prev and next point with the distances to the
	 * new point.
	 * 
	 * @param pts
	 * @param nw  new point that have to be inserted
	 * @param gf
	 * @return index (is then index of current point where cakeslice has to be
	 *         computed)
	 */
	public static int getIndexToInsertPoint(List<Point2D> pts, Point2D nw, GeometryFactory gf) {
		// find index where to add new point: minimize distances to neighbouring points
		// and compare them with straight distance -> must be equal
		int idx = 0;
		double minDist = Double.POSITIVE_INFINITY;
		// go through all consecutive points in ls (current pts (is added with each
		// inters point!))
		for (int o = 0; o < pts.size() - 1; o++) {
			Point prev = gf.createPoint(new Coordinate(pts.get(o).getX(), pts.get(o).getY()));
			Point next = gf.createPoint(new Coordinate(pts.get(o + 1).getX(), pts.get(o + 1).getY()));
			// distance between prev and next
			double dist_prev_next = prev.distance(next);
			// distances to new intersection point
			double dist_prev_new = prev.distance(gf.createPoint(new Coordinate(nw.getX(), nw.getY())));
			double dist_next_new = next.distance(gf.createPoint(new Coordinate(nw.getX(), nw.getY())));
			double dist_sum = dist_prev_new + dist_next_new;

			double tol = 0.000001;

			if (Math.abs(dist_prev_next - dist_sum) <= tol) {

//				System.out.println("distance prev_next: " + dist_prev_next);
//				System.out.println("distance sum prev-new and new-next: " + dist_sum);

				if (dist_sum <= minDist) {
					minDist = dist_sum;
					// store indices
					idx = o + 1; // was index of next, is now index of current
				}
			}
		}
		return idx;
	}

	/**
	 * this function creates a cakeslice of three parameters
	 * 
	 * @param r1  large radius
	 * @param r2  smaller radius (if 0 -> normal large cakesice)
	 * @param phi openig angle of cakeslice
	 */
	public static Polygon createCakeslice(double r1, double r2, double phi) {
		if (r2 != 0.0) {
			GeometryFactory gf = new GeometryFactory();
			GeometricShapeFactory shapeFactory = new GeometricShapeFactory(gf);
			shapeFactory.setCentre(new Coordinate(0, 0));
			// create coords of outer circle
			shapeFactory.setSize(2 * r1);
			LineString arc = shapeFactory.createArc(0, phi);
			Coordinate[] coords = new Coordinate[arc.getCoordinates().length + 2];
			coords[0] = new Coordinate(0, 0);
			for (int i = 0; i < arc.getCoordinates().length; ++i) {
				coords[i + 1] = arc.getCoordinateN(i);
			}
			coords[coords.length - 1] = new Coordinate(0, 0);
			// create coords of inner circle
			GeometryFactory gf2 = new GeometryFactory();
			GeometricShapeFactory shapeFactory2 = new GeometricShapeFactory(gf2);
			shapeFactory2.setSize(2 * r2);
			LineString arc2 = shapeFactory2.createArc(0, phi);
			Coordinate[] coords2 = new Coordinate[arc2.getCoordinates().length + 1];
			coords2[0] = new Coordinate(0, 0);
			for (int i = 0; i < arc2.getCoordinates().length; ++i) {
				coords2[i + 1] = arc2.getCoordinateN(i);
			}
			coords2[coords2.length - 1] = new Coordinate(0, 0);
			// merge all coords into one array and repeat first coordinate
			Coordinate[] coords_all = new Coordinate[arc.getCoordinates().length + arc2.getCoordinates().length + 1];
			for (int j = 0; j < coords.length; j++) {
				coords_all[j] = coords[j];
			}
			int o = 0;
			for (int j = coords.length; j < coords2.length; j++) {
				coords_all[j] = coords2[o];
				o++;
			}
			coords_all[coords_all.length - 1] = new Coordinate(0, 0);
			Polygon poly = gf.createPolygon(coords);
			return poly;
		} else {
			System.out.println("r2 = 0 -> normal cakeslice");
			GeometryFactory gf = new GeometryFactory();
			GeometricShapeFactory shapeFactory = new GeometricShapeFactory(gf);
			shapeFactory.setCentre(new Coordinate(0, 0));
			shapeFactory.setSize(2 * r1);
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

	public static final Polygon createEnvelope(double[] bounds) {
		return createEnvelope(bounds[0], bounds[1], bounds[2], bounds[3]);
	}

	public static final Polygon createEnvelope(double xMin, double yMin, double xMax, double yMax) {
		return createEnvelope(xMin, yMin, xMax, yMax, "EPSG:3857", "EPSG:25832");
	}

	public static final Polygon createEnvelope(double xMin, double yMin, double xMax, double yMax, String sourceCode,
			String targetCode) {
		Coordinate bottomLeft = new Coordinate(xMin, yMin);
		Coordinate bottomRight = new Coordinate(xMax, yMin);
		Coordinate topLeft = new Coordinate(xMin, yMax);
		Coordinate topRight = new Coordinate(xMax, yMax);

		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		Polygon polyInSourceCRS = gf
				.createPolygon(new Coordinate[] { bottomLeft, topLeft, topRight, bottomRight, bottomLeft });

		try {
			CoordinateReferenceSystem sourceCRS = CRS.decode(sourceCode);
			CoordinateReferenceSystem targetCRS = CRS.decode(targetCode);

			MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS);
			return (Polygon) JTS.transform(polyInSourceCRS, transform);
		} catch (NoSuchAuthorityCodeException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error in createEnvelope");
		} catch (FactoryException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error in createEnvelope");
		} catch (MismatchedDimensionException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error in createEnvelope");
		} catch (TransformException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error in createEnvelope");
		}
	}

}
