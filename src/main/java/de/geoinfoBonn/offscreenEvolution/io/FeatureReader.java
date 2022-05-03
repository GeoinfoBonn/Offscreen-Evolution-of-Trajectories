package de.geoinfoBonn.offscreenEvolution.io;

import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.net.MalformedURLException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.Property;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

/**
 * Reader to import all Features from an ESRI shapefile.
 * 
 * @author Axel Forsch
 *
 */
public class FeatureReader {
	private static final Logger logger = Logger.getLogger(MethodHandles.lookup().lookupClass().getName());

	/**
	 * Reads an ESRI Shapefile and returns a list of raw JTS-geometries without the
	 * features' attributes.
	 * 
	 * @param shapefile location of the shapefile (.shp) to read
	 * @return LinkedList of the geometries
	 * @throws IOException
	 */
	public static List<Geometry> readRawGeometryFromShapefile(File shapefile) throws IOException {
		long time = System.currentTimeMillis();

		Map<String, Object> map = new HashMap<>();
		map.put("url", shapefile.toURI().toURL());

		DataStore dataStore = DataStoreFinder.getDataStore(map);
		String typeName = dataStore.getTypeNames()[0];

		FeatureSource<SimpleFeatureType, SimpleFeature> source = dataStore.getFeatureSource(typeName);

		FeatureCollection<SimpleFeatureType, SimpleFeature> collection = source.getFeatures();

		LinkedList<Geometry> result = new LinkedList<>();

		try (FeatureIterator<SimpleFeature> features = collection.features()) {
			while (features.hasNext()) {
				SimpleFeature feature = features.next();

				Object geometry = feature.getDefaultGeometry();
				if (geometry instanceof Geometry) {
					result.add((Geometry) geometry);
				}
			}
		} finally {
			dataStore.dispose();
		}

		time = System.currentTimeMillis() - time;
		logger.fine("Finished reading shapefile, found " + result.size() + " geometries.");

		return result;
	}

	/**
	 * Reads an ESRI Shapefile and returns a list of {@link Feature} containing the
	 * geometry as well as all the attributes.
	 * 
	 * @param shapefile location of the shapefile (.shp) to read
	 * @return LinkedList of the features
	 */
	public static List<Feature> readFeatureFromShapefile(File shapefile) {
		return readFeatureFromShapefile(shapefile, null);
	}

	/**
	 * Reads an ESRI Shapefile and returns a list of {@link Feature} containing the
	 * geometry as well as all the attributes.
	 * 
	 * @param shapefile location of the shapefile (.shp) to read
	 * @return LinkedList of the features
	 */
	public static List<Feature> readFeatureFromShapefile(String shapefile) {
		return readFeatureFromShapefile(new File(shapefile), null);
	}

	/**
	 * Reads an ESRI Shapefile and returns a list of {@link Feature} containing the
	 * geometry as well as all the attributes. Additionally, an
	 * {@link FeatureProcessor} can be provided that performs an action on each read
	 * feature.
	 * 
	 * Here is an example for using this method to import a graph with data
	 * downloaded from Geofabrik using
	 * {@link de.geoinfoBonn.graphLibrary.io.shp.GeofabrikRoadGraphFactory
	 * GeofabrikRoadGraphFactory} and
	 * {@link de.geoinfoBonn.graphLibrary.io.shp.GraphProcessor GraphProcessor}
	 * 
	 * <pre>
	 * GraphFactory<BasicGeometricGraph<Point2D, DoubleWeightDataWithInfo<GeofabrikRoadData>>, Point2D, DoubleWeightDataWithInfo<GeofabrikRoadData>> factory = new GeofabrikRoadGraphFactory();
	 * GraphProcessor<BasicGeometricGraph<Point2D, DoubleWeightDataWithInfo<GeofabrikRoadData>>, Point2D, DoubleWeightDataWithInfo<GeofabrikRoadData>> processor = new GraphProcessor<>(
	 * 		factory);
	 *
	 * BasicGeometricGraph<Point2D, DoubleWeightDataWithInfo<GeofabrikRoadData>> graph;
	 * try {
	 * 	FeatureReader.readFeatureFromShapefile(new File(filename), processor);
	 * 	graph = processor.getGraph();	 *
	 * } catch (IOException e) {
	 * 	e.printStackTrace();
	 * }
	 * </pre>
	 * 
	 * @param shapefile location of the shapefile (.shp) to read
	 * @param processor defines an action to perform on each read feature
	 * @return LinkedList of the features
	 * @see <a href="download.geofabrik.de/">download.geofabrik.de/</a>
	 */
	public static List<Feature> readFeatureFromShapefile(File shapefile, FeatureProcessor processor) {
		long time = System.currentTimeMillis();
		try {

			Map<String, Object> map = new HashMap<>();

			map.put("url", shapefile.toURI().toURL());

			DataStore dataStore = DataStoreFinder.getDataStore(map);
			String typeName = dataStore.getTypeNames()[0];

			FeatureSource<SimpleFeatureType, SimpleFeature> source = dataStore.getFeatureSource(typeName);
			if (source.getSchema().getCoordinateReferenceSystem() != null)
				logger.finer("CRS of Shapefile: " + source.getSchema().getCoordinateReferenceSystem().getName());

			FeatureCollection<SimpleFeatureType, SimpleFeature> collection = source.getFeatures();
			Envelope extend = extendFromFeatureCollection(collection);
			if (processor != null)
				processor.initialize(extend, source.getSchema().getCoordinateReferenceSystem());

			LinkedList<Feature> result = new LinkedList<>();
			Feature currentFeature = null;
			String name;
			Object value;

			int perc = collection.size() / 100;
			int k = 0;
			int counter = 0;
			boolean verbose = false;

			try (FeatureIterator<SimpleFeature> features = collection.features()) {
				while (features.hasNext()) {
					if (verbose) {
						if (counter >= k * perc) {
							logger.finer(k + "% of features imported");
							k += 10;
						}
						++counter;
					}

					SimpleFeature feature = features.next();

					Object geometry = feature.getDefaultGeometry();
					try {
						currentFeature = new Feature((Geometry) geometry);
					} catch (ClassCastException e) {
						logger.warning("feature geometry has wrong type!");
						continue;
					}

					for (Property property : feature.getProperties()) {
						name = property.getName().getLocalPart();
						value = property.getValue();
						if (name != "the_geom") // geometry is already added separately, should not appear in attribute
												// list
							currentFeature.setAttribute(name, value);
					}

					if (processor != null)
						currentFeature = processor.processFeature(currentFeature);
					result.add(currentFeature);

				}
			} finally {
				dataStore.dispose();
			}

			time = System.currentTimeMillis() - time;
			logger.fine("Finished reading shapefile, found " + result.size() + " geometries.");
			return result;
		} catch (MalformedURLException e1) {
			e1.printStackTrace();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return null;
	}

	private static Envelope extendFromFeatureCollection(FeatureCollection<SimpleFeatureType, SimpleFeature> fc) {
		ReferencedEnvelope re = fc.getBounds();
		return new Envelope(re.getMinX(), re.getMaxX(), re.getMinY(), re.getMaxY());
	}
}
