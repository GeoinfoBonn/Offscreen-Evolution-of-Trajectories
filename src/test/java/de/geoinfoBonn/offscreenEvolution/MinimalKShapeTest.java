package de.geoinfoBonn.offscreenEvolution;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.stream.Stream;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import de.geoinfoBonn.offscreenEvolution.shapefinder.MinimalKShapeFinder;
import de.geoinfoBonn.offscreenEvolution.shapefinder.RectanglePreparer;

/**
 * Unit test for simple App.
 */
public class MinimalKShapeTest {
	private ArrayList<Point2D> points;

	@BeforeEach
	public void initializePoints() {
		points = new ArrayList<>();
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
	}

	public static Stream<Arguments> rectangleProvider() {
		return Stream.of(//
				Arguments.of(1, true, 2.0, 2.0), //
				Arguments.of(1, true, 2.0, 2.0), //
				Arguments.of(2, true, 4.0, 6.0), //
				Arguments.of(2, true, 4.0, 6.0), //
				Arguments.of(3, true, 7.0, 6.0), //
				Arguments.of(3, true, 7.0, 6.0), //
				Arguments.of(4, true, 7.0, 12.0), //
				Arguments.of(4, false, 4.0, 17.0), //
				Arguments.of(5, true, 14.0, 8.0), //
				Arguments.of(5, false, 14.0, 8.0), //
				Arguments.of(6, true, 14.0, 11.0), //
				Arguments.of(6, false, 19.0, 8.0), //
				Arguments.of(7, true, 14.0, 12.0), //
				Arguments.of(7, false, 14.0, 12.0), //
				Arguments.of(8, true, 19.0, 12.0), //
				Arguments.of(8, false, 19.0, 12.0), //
				Arguments.of(9, true, 14.0, 18.0), //
				Arguments.of(9, false, 21.0, 12.0), //
				Arguments.of(10, true, 24.0, 12.0), //
				Arguments.of(10, false, 24.0, 12.0) //
		);
	}

	@ParameterizedTest
	@MethodSource("rectangleProvider")
	public void minRectangle(int k, boolean usePerimeter, double x, double y) {
		double[] res = MinimalKShapeFinder.findMinimalKShape(points, k, new RectanglePreparer(usePerimeter));
		Assertions.assertEquals(x, res[0]);
		Assertions.assertEquals(y, res[1]);
	}

}
