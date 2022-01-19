package de.geoinfoBonn.offscreenEvolution.shapefinder;

import java.util.Objects;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.util.GeometricShapeFactory;

public class ShapeCreator {

	/**
	 * Creates a geometry representing staggered cake slice sections. For each of
	 * the {@code n} parameters in {@code gamma}/{@code phi}/{@code r} one section
	 * gets created. Returns a geometry of type Polygon or MultiPolygon without
	 * holes.
	 * 
	 * @param centerPoint center point of the cake slice formation
	 * @param gamma       angle between x-axis and right border of the i-th cake
	 *                    slice, one value per section
	 * @param phi         opening angle of the i-th cake slice
	 * @param r           radius of the outer border of the i-th cake slice,
	 *                    measured from the center point, one value per section
	 * @return Polygon or MultiPolygon representing the staggered cake slice
	 *         formation, all geometries are without holes
	 */
	public static Geometry createStaggeredCakeslices(Coordinate centerPoint, double[] gamma, double[] phi, double[] r) {
		Objects.requireNonNull(centerPoint);
		Objects.requireNonNull(gamma);
		Objects.requireNonNull(phi);
		Objects.requireNonNull(r);

		// number of cakeslice (and cakesliceSections) to draw
		int n = gamma.length;

		// check input
		if (phi.length != n || r.length != n)
			throw new IllegalArgumentException("length of gamma, phi and r must be the same! Found: #gamma="
					+ gamma.length + ", #phi=" + phi.length + ", #r=" + r.length);

		Polygon[] polygons = new Polygon[n];

		double r1 = 0;
		for (int i = 0; i < n; r1 = r[i], ++i) {
			Polygon currentCakeslice = createCakesliceSection(centerPoint, gamma[i], phi[i], r[i], r1);
			polygons[i] = currentCakeslice;
		}

		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		GeometryCollection collection = gf.createGeometryCollection(polygons);

		// buffering has two effects:
		// 1) polygons get unioned
		// 2) small holes between polygons which occur due to numeric issues in the arc
		// creation get closed.
		BufferParameters bufferParam = new BufferParameters(//
				BufferParameters.DEFAULT_QUADRANT_SEGMENTS, //
				BufferParameters.CAP_FLAT, //
				BufferParameters.JOIN_MITRE, //
				BufferParameters.DEFAULT_MITRE_LIMIT);

		double bufferSize = 1e-1;
		Geometry result = BufferOp.bufferOp(collection, bufferSize, bufferParam);
		result = BufferOp.bufferOp(result, -bufferSize, bufferParam);

		// ##### check if assumptions on output are correct

		// first assumption: result geometry is of type Polygon or MultiPolygon
		if (!result.getGeometryType().endsWith("Polygon"))
			System.err.println(
					"The returned staggered cakeslice is neither a Polygon nor an MultiPolygon, please investigate why this is the case!");

		// second assumption: result geometries have no holes
		for (int i = 0; i < result.getNumGeometries(); i++) {
			int numHoles = ((Polygon) result.getGeometryN(i)).getNumInteriorRing();
			if (numHoles > 0)
				System.err.println("Geometry " + i + " of the result has " + numHoles
						+ " holes, please investigate why this is the case! Maybe the buffer size needs to be enlarged.");
		}

		return result;
	}

	/**
	 * Creates a geometry representing staggered cake slice sections. For each of
	 * the {@code n} parameters in {@code gamma}/{@code phi}/{@code r} one section
	 * gets created. Returns a geometry of type Polygon or MultiPolygon without
	 * holes. The Geometry is scaled with factor s.
	 * 
	 * @param centerPoint center point of the cake slice formation
	 * @param gamma       angle between x-axis and right border of the i-th cake
	 *                    slice, one value per section
	 * @param phi         opening angle of the i-th cake slice
	 * @param r           radius of the outer border of the i-th cake slice,
	 *                    measured from the center point, one value per section
	 * @param s           scaling factor
	 * @return Scaled Polygon or MultiPolygon representing the staggered cake slice
	 *         formation, all geometries are without holes
	 */
	public static Geometry createStaggeredCakeslices(Coordinate centerPoint, double[] gamma, double[] phi, double[] r,
			double s) {
		Objects.requireNonNull(centerPoint);
		Objects.requireNonNull(gamma);
		Objects.requireNonNull(phi);
		Objects.requireNonNull(r);

		double[] r_scaled = new double[r.length];
		for (int i = 0; i < r.length; i++) {
			r_scaled[i] = r[i] * s;
		}

		return createStaggeredCakeslices(centerPoint, gamma, phi, r_scaled);
	}

	/**
	 * Creates a cake slice section with the given parameters and returns a simple
	 * Polygon.
	 * 
	 * @param centerPoint center point of the cake slice
	 * @param gamma       angle between x-axis and right border of the cake slice
	 * @param phi         opening angle of the cake slice
	 * @param outerRadius radius of the outer border of the cake slice, measured
	 *                    from the center point, must be non-negative
	 * @param innerRadius radius of the inner border of the cake slice, measured
	 *                    from the center point, must be non-negative, value of 0
	 *                    means a whole cake slice instead of just a section gets
	 *                    created
	 * @return a simple Polygon representing the cake slice
	 */
	public static Polygon createCakesliceSection(Coordinate centerPoint, double gamma, double phi, double outerRadius,
			double innerRadius) {
		Objects.requireNonNull(centerPoint);
		if (outerRadius < 0 || innerRadius < 0)
			throw new IllegalArgumentException(
					"Radii must be positive, found: outerRadius=" + outerRadius + ", innerRadius=" + innerRadius);

		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();

		GeometricShapeFactory shapeFactory = new GeometricShapeFactory(gf);
		shapeFactory.setCentre(centerPoint);
		Geometry innerArc = null;
		if (innerRadius != 0) {
			shapeFactory.setSize(innerRadius);
			innerArc = shapeFactory.createArc(gamma, phi);
		} else {
			innerArc = gf.createPoint(centerPoint);
		}
		shapeFactory.setSize(outerRadius);
		LineString outerArc = shapeFactory.createArc(gamma, phi).reverse();

		// outer shell of cakeslice, +1 as first point needs to be repeated to close
		// polygon
		Coordinate[] shell = new Coordinate[innerArc.getNumPoints() + outerArc.getNumPoints() + 1];
		for (int i = 0; i < innerArc.getNumPoints(); i++) {
			shell[i] = innerArc.getCoordinates()[i];
		}
		for (int i = 0, j = innerArc.getNumPoints(); i < outerArc.getNumPoints(); ++i, ++j) {
			shell[j] = outerArc.getCoordinateN(i);
		}
		shell[shell.length - 1] = innerArc.getCoordinates()[0]; // repeat first coordinate

		// create Polygon
		return gf.createPolygon(shell);
	}

	/**
	 * Creates a cake slice with the given parameters and returns a simple Polygon.
	 * 
	 * @param centerPoint center point of the cake slice
	 * @param gamma       angle between x-axis and right border of the cake slice
	 * @param phi         opening angle of the cake slice
	 * @param r           radius of the outer border of the cake slice, measured
	 *                    from the center point
	 * @return
	 */
	public static Polygon createCakeslice(Coordinate centerPoint, double gamma, double phi, double r) {
		Objects.requireNonNull(centerPoint);
		return createCakesliceSection(centerPoint, gamma, phi, r, 0);
	}

	/**
	 * Creates a scaled cake slice with the given parameters and returns a simple
	 * Polygon.
	 * 
	 * @param centerPoint center point of the cake slice
	 * @param gamma       angle between x-axis and right border of the cake slice
	 * @param phi         opening angle of the cake slice
	 * @param r           radius of the outer border of the cake slice, measured
	 *                    from the center point
	 * @param s           scaling factor
	 * @return
	 */
	public static Polygon createCakeslice(Coordinate centerPoint, double gamma, double phi, double r, double s) {
		Objects.requireNonNull(centerPoint);
		double r_scaled = r * s;
		return createCakesliceSection(centerPoint, gamma, phi, r_scaled, 0);
	}
}
