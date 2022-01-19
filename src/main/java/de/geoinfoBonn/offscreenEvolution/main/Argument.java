package de.geoinfoBonn.offscreenEvolution.main;

import java.lang.reflect.Field;
import java.util.Iterator;
import java.util.logging.Logger;

public class Argument {
	private static final Logger LOGGER = Logger.getLogger(Argument.class.getName());

	public static final String HELP = "-h";
	public static final String TRACK = "-in";
	public static final String FRAME_FILE = "-frame";
	public static final String N = "-n";
	public static final String K = "-k";
	public static final String KAPPA = "-kappa";
	public static final String B = "-b";
	public static final String BETA = "-beta";
	public static final String NUM_SHAPES = "-ns";
	public static final String TIME = "-t";
	public static final String DISTANCE = "-d";
	public static final String RANGE = "-r";
	public static final String OVERLAP = "-ol";
	public static final String VARIANT = "-var";
	public static final String SCALE_FACTOR = "-s";
	public static final String MAX_SPEED = "-vmax";
	public static final String OUTPUT = "-out";
	public static final String FRAME_BOUNDS = "-bounds";

	public static String getHelpText(String argumentIdentifier) {
		switch (argumentIdentifier) {
		case HELP:
			return argumentIdentifier + "\t" + "Prints a help text for this app.";
		case TRACK:
			return argumentIdentifier + "\t" + "Path to the trajectory shapefile (given in EPSG:xxx).";
		case FRAME_FILE:
			return argumentIdentifier + "\t" + "Path to the frame shapefile (given in EPSG:xxx).";
		case N:
			return argumentIdentifier + "\t" + "Number of points to look ahead [-]. (NOT IMPLEMENTED YET!)";
		case K:
			return argumentIdentifier + "\t" + "Number of points inside shape [-]. (NOT IMPLEMENTED YET!)";
		case KAPPA:
			return argumentIdentifier + "\t" + "Share of points inside shape [rel] (overrides '" + K + "').";
		case B:
			return argumentIdentifier + "\t" + "Width of border area [m].";
		case BETA:
			return argumentIdentifier + "\t" + "Share of map used for border [rel] (overrides '" + B + "').";
		case NUM_SHAPES:
			return argumentIdentifier + "\t" + "Number of shapes per signature [-].";
		case TIME:
			return argumentIdentifier + "\t" + "Time considered for each shape of the signature [min].";
		case OVERLAP:
			return argumentIdentifier + "\t" + "Temporal overlap between successive shapes of a signature [rel].";
		case VARIANT:
			return argumentIdentifier + "\t" + "Choose variant for signatures.";
		case SCALE_FACTOR:
			return argumentIdentifier + "\t" + "Scale factor of the signatures [rel].";
		case MAX_SPEED:
			return argumentIdentifier + "\t" + "Maximum speed of trajectories [m/s] (overrides '" + SCALE_FACTOR
					+ "').";
		case OUTPUT:
			return argumentIdentifier + "\t" + "Directory to write results to.";
		default:
			return argumentIdentifier + "\t" + "Identifier not found!";
		}
	}

	public static Iterator<String> iterator() {
		return new Iterator<String>() {
			Field[] fields = Argument.class.getFields();
			int i = 0;

			@Override
			public boolean hasNext() {
				return i < fields.length;
			}

			@Override
			public String next() {
				String ret;
				try {
					ret = (String) fields[i].get(null);
				} catch (IllegalArgumentException | IllegalAccessException e) {
					e.printStackTrace();
					return null;
				}
				++i;
				return ret;
			}
		};
	}

	public static String getOptionalArg(String[] args, String identifier) {
		for (int i = 0; i < args.length - 1; i++) {
			if (args[i].equals(identifier)) {
				LOGGER.finest("Read argument #" + (i + 1) + ": " + args[i + 1]);
				return args[i + 1];
			}
		}
		String errorMsg = "Mandatory argument '" + identifier + "' not found.";
		LOGGER.severe(errorMsg);
		throw new RuntimeException(errorMsg);
	}

	public static boolean containsOptionalArg(String[] args, String identifier) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals(identifier))
				return true;
		}
		return false;
	}
}