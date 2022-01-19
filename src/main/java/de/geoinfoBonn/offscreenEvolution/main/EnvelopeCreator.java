package de.geoinfoBonn.offscreenEvolution.main;

import java.io.File;
import java.util.LinkedList;

import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import de.geoinfoBonn.offscreenEvolution.utils.Tools;

public class EnvelopeCreator {

	public static void main(String[] args) {
		// coordinates given in crsFrom
		double xMin = Double.parseDouble(args[0]);
		double yMin = Double.parseDouble(args[1]);
		double xMax = Double.parseDouble(args[2]);
		double yMax = Double.parseDouble(args[3]);

		// output files
		File outputPath = new File(args[4]);

		// read codes of crs
		String targetCode = args.length > 5 ? args[5] : "EPSG:25832";
		String sourceCode = args.length > 6 ? args[6] : "EPSG:3857";

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
			Polygon polyInTargetCRS = (Polygon) JTS.transform(polyInSourceCRS, transform);

			LinkedList<Geometry> source = new LinkedList<>();
			source.add(polyInSourceCRS);

			LinkedList<Geometry> target = new LinkedList<>();
			target.add(polyInTargetCRS);

			Tools.exportGeometryAsPolygons(new File(outputPath, "envelopeInSource.shp").getAbsolutePath(), source, null,
					null, sourceCRS);
			Tools.exportGeometryAsPolygons(new File(outputPath, "envelopeInTarget.shp").getAbsolutePath(), target, null,
					null, targetCRS);

		} catch (NoSuchAuthorityCodeException e) {
			e.printStackTrace();
		} catch (FactoryException e) {
			e.printStackTrace();
		} catch (MismatchedDimensionException e) {
			e.printStackTrace();
		} catch (TransformException e) {
			e.printStackTrace();
		}

	}
}
