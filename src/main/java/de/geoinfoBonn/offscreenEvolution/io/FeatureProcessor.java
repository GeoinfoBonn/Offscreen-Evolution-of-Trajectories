package de.geoinfoBonn.offscreenEvolution.io;

import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Feature processors can be used to create derived data structures from
 * Shapefiles during the read process. The processors
 * {@link processFeature(Feature)} method is called on each Feature read by
 * {@link de.geoinfoBonn.graphLibrary.io.shp.FeatureReader#readFeatureFromShapefile(File, FeatureProcessor)}.
 * 
 * The method {@link initialize(Envelope)} is called once before reading the
 * Features and gets the Features extent as parameter. This can for example be
 * used to initialize QuadTrees for all Features.
 * 
 * @author Axel Forsch
 * 
 */
public interface FeatureProcessor {
	/**
	 * Called once before reading the Features, this method can be used to
	 * initialize anything needing the extent of all features.
	 * 
	 * @param env extent of all Features in the file
	 */
	public void initialize(Envelope env, CoordinateReferenceSystem coordinateReferenceSystem);

	/**
	 * This method is called on each Feature in the input file. This method can be
	 * used to build a derived data structure from the Features (such as a e.g. a
	 * graph). The Feature can be altered inside this method and retured to modify
	 * the list of Features
	 * {@link de.geoinfoBonn.graphLibrary.io.shp.FeatureReader#readFeatureFromShapefile(File, FeatureProcessor)}
	 * returns.
	 * 
	 * @param f read Feature
	 * @return (possibly modified) Feature to put into result of
	 *         {@link de.geoinfoBonn.graphLibrary.io.shp.FeatureReader#readFeatureFromShapefile(File, FeatureProcessor)}
	 */
	public Feature processFeature(Feature f);
}