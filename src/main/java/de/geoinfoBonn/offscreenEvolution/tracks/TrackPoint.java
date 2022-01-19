package de.geoinfoBonn.offscreenEvolution.tracks;

import java.time.LocalDateTime;

import org.locationtech.jts.geom.Coordinate;

public class TrackPoint implements Comparable<TrackPoint> {

	public enum PointType {
		TRACK_POINT, // normal track point
		CROSSING_FORWARD, // crossing point where the outside is forward
		CROSSING_BACKWARD, // crossing point where the outside is backward
		CROSSING_OUTER
	}

	private Coordinate location;
	private LocalDateTime timestamp;
	private double sequence;
	private PointType type;

	public TrackPoint(Coordinate location, LocalDateTime timestamp, double sequence) {
		this.location = location;
		this.timestamp = timestamp;
		this.sequence = sequence;
		this.type = PointType.TRACK_POINT;
	}

	public TrackPoint(Coordinate location, LocalDateTime timestamp, double sequence, PointType type) {
		this.location = location;
		this.timestamp = timestamp;
		this.sequence = sequence;
		this.type = type;
	}

	public double getSequence() {
		return sequence;
	}

	public double getSequence(double minSequence) {
		return sequence - minSequence;
	}

	public LocalDateTime getTimestamp() {
		return timestamp;
	}

	public Coordinate getCoordinate() {
		return location;
	}

	@Override
	public int compareTo(TrackPoint o) {
		return (int) Math.signum(sequence - o.sequence);
	}

	public PointType getType() {
		return type;
	}

	public void setType(PointType type) {
		this.type = type;
	}

	@Override
	public String toString() {
		return type + " [location=" + location + ", timestamp=" + timestamp + ", sequence=" + sequence + "]";
	}
}
