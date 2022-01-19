package de.geoinfoBonn.offscreenEvolution.main;

import java.io.File;
import java.util.Optional;
import java.util.logging.Logger;

public class RunConfig {
	private static final Logger LOGGER = Logger.getLogger(RunConfig.class.getName());

	private String idColumn;
	private String segmentColumn;
	private String seqColumn;
	private String timeColumn;
	private String timePattern;

	private File trackFile;
	private Optional<File> frameFile;
	private Optional<double[]> pseudomercatorFrameBounds;
	private File outputDir;
	private Optional<Integer> b;
	private Optional<Double> beta;
	private Optional<Integer> k;
	private Optional<Double> kappa;
	private Optional<Integer> n;
	private Optional<Integer> maxTime;
	private Optional<Integer> maxDist;
	private Optional<Integer> maxRange;
	private int numShapes;
	private double temporalOverlap;
	private Optional<Double> scaleFactor;
	private Optional<Double> maxSpeed;
	private Variant variant;

	private RunConfig(String idColumn, String segmentColumn, String seqColumn, String timeColumn, String timePattern,
			File trackFile, Optional<File> frameFile, Optional<double[]> pseudomercatorFrameBounds, File outputDir,
			Optional<Integer> b, Optional<Double> beta, Optional<Integer> k, Optional<Double> kappa,
			Optional<Integer> n, Optional<Integer> maxTime, Optional<Integer> maxDist, Optional<Integer> maxRange,
			int numShapes, double temporalOverlap, Optional<Double> scaleFactor, Optional<Double> maxSpeed,
			Variant variant) {
		super();
		this.idColumn = idColumn;
		this.segmentColumn = segmentColumn;
		this.seqColumn = seqColumn;
		this.timeColumn = timeColumn;
		this.timePattern = timePattern;
		this.trackFile = trackFile;
		this.frameFile = frameFile;
		this.pseudomercatorFrameBounds = pseudomercatorFrameBounds;
		this.outputDir = outputDir;
		this.b = b;
		this.beta = beta;
		this.k = k;
		this.kappa = kappa;
		this.n = n;
		this.maxTime = maxTime;
		this.maxDist = maxDist;
		this.maxRange = maxRange;
		this.numShapes = numShapes;
		this.temporalOverlap = temporalOverlap;
		this.scaleFactor = scaleFactor;
		this.maxSpeed = maxSpeed;
		this.variant = variant;
	}

	public void logRunConfig() {
		LOGGER.config("idColum: " + idColumn);
		LOGGER.config("segmentColumn: " + segmentColumn);
		LOGGER.config("seqColumn: " + seqColumn);
		LOGGER.config("timeColumn: " + timeColumn);
		LOGGER.config("timePattern: " + timePattern);
		LOGGER.config("trackFile: " + trackFile);
		LOGGER.config("frameFile: " + frameFile);
		LOGGER.config("b: " + (b.isPresent() ? b.get() : "not set"));
		LOGGER.config("beta: " + (beta.isPresent() ? beta.get() : "not set"));
		LOGGER.config("k: " + (k.isPresent() ? k.get() : "not set"));
		LOGGER.config("kappa: " + (kappa.isPresent() ? kappa.get() : "not set"));
		LOGGER.config("n: " + (n.isPresent() ? n.get() : "not set"));
		LOGGER.config("maxTime: " + (maxTime.isPresent() ? maxTime.get() : "not set"));
		LOGGER.config("maxDist: " + (maxDist.isPresent() ? maxDist.get() : "not set"));
		LOGGER.config("maxRange: " + (maxRange.isPresent() ? maxRange.get() : "not set"));
		LOGGER.config("numShapes: " + numShapes);
		LOGGER.config("temporalOverlap: " + temporalOverlap);
		LOGGER.config("scaleFactor: " + (scaleFactor.isPresent() ? scaleFactor.get() : "not set"));
		LOGGER.config("maxSpeed: " + (maxSpeed.isPresent() ? maxSpeed.get() : "not set"));
		LOGGER.config("variant: " + variant);
	}

	public static RunConfig parseArguments(String[] args) {
		RunConfigBuilder builder = RunConfig.builder();

		if (Argument.containsOptionalArg(args, Argument.BETA)) {
			double beta = Double.parseDouble(Argument.getOptionalArg(args, Argument.BETA));
			builder.beta(beta);
		} else if (Argument.containsOptionalArg(args, Argument.B)) {
			int b = Integer.parseInt(Argument.getOptionalArg(args, Argument.B));
			builder.b(b);
		} else {
			String errorMsg = "Either '" + Argument.B + "' or '" + Argument.BETA + "' must be set.";
			LOGGER.severe(errorMsg);
			throw new RuntimeException(errorMsg);
		}

		if (Argument.containsOptionalArg(args, Argument.TIME)) {
			int maxTime = Integer.parseInt(Argument.getOptionalArg(args, Argument.TIME));
			builder.maxTime(maxTime);
		} else if (Argument.containsOptionalArg(args, Argument.N)) {
			int n = Integer.parseInt(Argument.getOptionalArg(args, Argument.N));
			builder.n(n);
		} else if (Argument.containsOptionalArg(args, Argument.DISTANCE)) {
			int maxDist = Integer.parseInt(Argument.getOptionalArg(args, Argument.DISTANCE));
			builder.maxDist(maxDist);
		} else if (Argument.containsOptionalArg(args, Argument.RANGE)) {
			int maxRange = Integer.parseInt(Argument.getOptionalArg(args, Argument.RANGE));
			builder.maxRange(maxRange);
		} else {
			String errorMsg = "Either '" + Argument.N + "' or '" + Argument.TIME + "' must be set.";
			LOGGER.severe(errorMsg);
			throw new RuntimeException(errorMsg);
		}

		if (Argument.containsOptionalArg(args, Argument.KAPPA)) {
			double kappa = Double.parseDouble(Argument.getOptionalArg(args, Argument.KAPPA));
			builder.kappa(kappa);
		} else if (Argument.containsOptionalArg(args, Argument.K)) {
			int k = Integer.parseInt(Argument.getOptionalArg(args, Argument.K));
			builder.k(k);
		} else {
			String errorMsg = "Either '" + Argument.K + "' or '" + Argument.KAPPA + "' must be set.";
			LOGGER.severe(errorMsg);
			throw new RuntimeException(errorMsg);
		}

		if (Argument.containsOptionalArg(args, Argument.NUM_SHAPES)) {
			int numShapes = Integer.parseInt(Argument.getOptionalArg(args, Argument.NUM_SHAPES));
			builder.numShapes(numShapes);
		}

		if (Argument.containsOptionalArg(args, Argument.OVERLAP)) {
			double temporalOverlap = Double.parseDouble(Argument.getOptionalArg(args, Argument.OVERLAP));
			builder.temporalOverlap(temporalOverlap);
		}

		if (Argument.containsOptionalArg(args, Argument.MAX_SPEED)) {
			double maxSpeed = Double.parseDouble(Argument.getOptionalArg(args, Argument.MAX_SPEED));
			builder.maxSpeed(maxSpeed);
		} else if (Argument.containsOptionalArg(args, Argument.SCALE_FACTOR)) {
			double scaleFactor = Double.parseDouble(Argument.getOptionalArg(args, Argument.SCALE_FACTOR));
			builder.scaleFactor(scaleFactor);
		} else if (Argument.containsOptionalArg(args, Argument.TIME)) {
			String errorMsg = "Either '" + Argument.SCALE_FACTOR + "' or '" + Argument.MAX_SPEED + "' must be set.";
			LOGGER.severe(errorMsg);
			throw new RuntimeException(errorMsg);
		}

		if (Argument.containsOptionalArg(args, Argument.VARIANT)) {
			Variant variant = Variant.valueOf(Argument.getOptionalArg(args, Argument.VARIANT).toUpperCase());
			builder.variant(variant);
		}

		if (Argument.containsOptionalArg(args, Argument.FRAME_FILE)) {
			String frameFilename = Argument.getOptionalArg(args, Argument.FRAME_FILE);
			builder.frameFile(frameFilename);
		} else if (Argument.containsOptionalArg(args, Argument.FRAME_BOUNDS)) {
			String bound = Argument.getOptionalArg(args, Argument.FRAME_BOUNDS);
			String[] bounds = bound.split(",");
			double[] boundss = new double[bounds.length];
			for (int i = 0; i < bounds.length; i++)
				boundss[i] = Double.parseDouble(bounds[i]);
			builder.pseudomercatorFrameBounds(boundss);
		}

		String trackFilename = Argument.getOptionalArg(args, Argument.TRACK).replaceAll("%20", " ");
		String outputDirname = Argument.getOptionalArg(args, Argument.OUTPUT).replaceAll("%20", " ");

		return builder.build(new File(trackFilename), new File(outputDirname));
	}

	public String getIdColumn() {
		return idColumn;
	}

	public String getSegmentColumn() {
		return segmentColumn;
	}

	public String getSeqColumn() {
		return seqColumn;
	}

	public String getTimeColumn() {
		return timeColumn;
	}

	public String getTimePattern() {
		return timePattern;
	}

	public File getTrackFile() {
		return trackFile;
	}

	public Optional<File> getFrameFile() {
		return frameFile;
	}

	public Optional<double[]> getPseudomercatorFrameBounds() {
		return pseudomercatorFrameBounds;
	}

	public File getOutputDir() {
		return outputDir;
	}

	public Optional<Integer> getB() {
		return b;
	}

	public Optional<Double> getBeta() {
		return beta;
	}

	public Optional<Integer> getK() {
		return k;
	}

	public Optional<Double> getKappa() {
		return kappa;
	}

	public Optional<Integer> getN() {
		return n;
	}

	public Optional<Integer> getMaxTime() {
		return maxTime;
	}

	public Optional<Integer> getMaxDist() {
		return maxDist;
	}

	public Optional<Integer> getMaxRange() {
		return maxRange;
	}

	public int getNumShapes() {
		return numShapes;
	}

	public double getOverlap() {
		return temporalOverlap;
	}

	public Optional<Double> getScaleFactor() {
		return scaleFactor;
	}

	public Optional<Double> getMaxSpeed() {
		return maxSpeed;
	}

	public Variant getVariant() {
		return variant;
	}

	public static RunConfigBuilder builder() {
		return new RunConfigBuilder();
	}

	public static class RunConfigBuilder {
		private String idColumn = "track_id";
		private String segmentColumn = "segment_id";
		private String seqColumn = "point_id";
		private String timeColumn = "time";
		private String timePattern = "yyyy/MM/dd HH:mm:ss.SSS";

		private Optional<File> frameFile;
		private Optional<double[]> pseudomercatorFrameBounds;

		private Optional<Integer> b;
		private Optional<Double> beta;
		private Optional<Integer> k;
		private Optional<Double> kappa;
		private Optional<Integer> n;
		private Optional<Integer> maxTime;
		private Optional<Integer> maxDist;
		private Optional<Integer> maxRange;
		private int numShapes = 1;
		private double temporalOverlap = 0;
		private Optional<Double> scaleFactor = Optional.empty();
		private Optional<Double> maxSpeed = Optional.empty();
		private Variant variant = Variant.ALL_POINTS;

		public RunConfig build(File trackFile, File outputDir) {
			if (frameFile.isEmpty() && pseudomercatorFrameBounds.isEmpty())
				throw new IllegalArgumentException("Either frameFile or pseudomercatorFrameBounds must be set");
			if (b.isEmpty() && beta.isEmpty())
				throw new IllegalArgumentException("Either b or beta must be set");
			if (k.isEmpty() && kappa.isEmpty())
				throw new IllegalArgumentException("Either k or kappa must be set");
			if (n.isEmpty() && maxTime.isEmpty() && maxDist.isEmpty() && maxRange.isEmpty())
				throw new IllegalArgumentException("Either n or maxTime or maxDist or maxRange must be set");
			if (maxTime.isPresent() && scaleFactor.isEmpty() && maxSpeed.isEmpty())
				throw new IllegalArgumentException("Either scaleFactor or maxSpeed must be set when t is used");
			if (k.isPresent() && n.isPresent() && k.get() > n.get())
				throw new IllegalArgumentException("k must be smaller or equal to n");
			if (maxSpeed.isPresent() && maxTime.isEmpty())
				throw new IllegalArgumentException("maxSpeed can only be used in conjunction with maxTime");

			return new RunConfig(idColumn, segmentColumn, seqColumn, timeColumn, timePattern, trackFile, frameFile,
					pseudomercatorFrameBounds, outputDir, b, beta, k, kappa, n, maxTime, maxDist, maxRange, numShapes,
					temporalOverlap, scaleFactor, maxSpeed, variant);
		}

		public RunConfigBuilder idColumn(String idColumn) {
			this.idColumn = idColumn;
			return this;
		}

		public RunConfigBuilder segmentColumn(String segmentColumn) {
			this.segmentColumn = segmentColumn;
			return this;
		}

		public RunConfigBuilder seqColumn(String seqColumn) {
			this.seqColumn = seqColumn;
			return this;
		}

		public RunConfigBuilder timeColumn(String timeColumn) {
			this.timeColumn = timeColumn;
			return this;
		}

		public RunConfigBuilder timePattern(String timePattern) {
			this.timePattern = timePattern;
			return this;
		}

		public RunConfigBuilder frameFile(String frameFileName) {
			this.frameFile = Optional.of(new File(frameFileName));
			this.pseudomercatorFrameBounds = Optional.empty();
			return this;
		}

		public RunConfigBuilder pseudomercatorFrameBounds(double[] bounds) {
			this.frameFile = Optional.empty();
			this.pseudomercatorFrameBounds = Optional.of(bounds);
			return this;
		}

		public RunConfigBuilder b(int b) {
			if (b <= 0) {
				String errorMsg = "b must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.b = Optional.of(b);
			this.beta = Optional.empty();
			return this;
		}

		public RunConfigBuilder beta(double beta) {
			if (beta <= 0 || beta >= 0.5) {
				String errorMsg = "beta must be in interval (0,0.5).";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.b = Optional.empty();
			this.beta = Optional.of(beta);
			return this;
		}

		public RunConfigBuilder k(int k) {
			if (k <= 0) {
				String errorMsg = "k must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.k = Optional.of(k);
			this.kappa = Optional.empty();
			return this;
		}

		public RunConfigBuilder kappa(double kappa) {
			if (kappa <= 0 || kappa > 1) {
				String errorMsg = "kappa must be in interval (0, 1].";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.k = Optional.empty();
			this.kappa = Optional.of(kappa);
			return this;
		}

		public RunConfigBuilder n(int n) {
			if (n <= 0) {
				String errorMsg = "n must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.n = Optional.of(n);
			this.maxTime = Optional.empty();
			this.maxDist = Optional.empty();
			return this;
		}

		public RunConfigBuilder maxTime(int maxTime) {
			if (maxTime <= 0) {
				String errorMsg = "Maximum time must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.n = Optional.empty();
			this.maxTime = Optional.of(maxTime);
			this.maxDist = Optional.empty();
			this.maxRange = Optional.empty();
			return this;
		}

		public RunConfigBuilder maxDist(int maxDist) {
			if (maxDist <= 0) {
				String errorMsg = "Maximum distance must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.n = Optional.empty();
			this.maxTime = Optional.empty();
			this.maxDist = Optional.of(maxDist);
			this.maxRange = Optional.empty();
			return this;
		}

		public RunConfigBuilder maxRange(int maxRange) {
			if (maxRange <= 0) {
				String errorMsg = "Maximum range must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.n = Optional.empty();
			this.maxTime = Optional.empty();
			this.maxDist = Optional.empty();
			this.maxRange = Optional.of(maxRange);
			return this;
		}

		public RunConfigBuilder numShapes(int numShapes) {
			if (numShapes <= 0) {
				String errorMsg = "Number of shapes must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.numShapes = numShapes;
			return this;
		}

		public RunConfigBuilder temporalOverlap(double temporalOverlap) {
			if (temporalOverlap < 0 || temporalOverlap >= 1) {
				String errorMsg = "Temporal overlap must be in interal [0, 1).";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.temporalOverlap = temporalOverlap;
			return this;
		}

		public RunConfigBuilder scaleFactor(double scaleFactor) {
			if (scaleFactor <= 0 || scaleFactor > 1) {
				String errorMsg = "Scale factor must be in interal (0, 1].";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.scaleFactor = Optional.of(scaleFactor);
			this.maxSpeed = Optional.empty();
			return this;
		}

		public RunConfigBuilder maxSpeed(double maxSpeed) {
			if (maxSpeed <= 0) {
				String errorMsg = "Max speed must be positive.";
				LOGGER.severe(errorMsg);
				throw new RuntimeException(errorMsg);
			}
			this.scaleFactor = Optional.empty();
			this.maxSpeed = Optional.of(maxSpeed);
			return this;
		}

		public RunConfigBuilder variant(Variant variant) {
			this.variant = variant;
			return this;
		}
	}
}
