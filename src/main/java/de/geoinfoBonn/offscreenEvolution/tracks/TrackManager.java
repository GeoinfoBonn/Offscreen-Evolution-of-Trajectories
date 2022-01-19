package de.geoinfoBonn.offscreenEvolution.tracks;

import java.util.ArrayList;
import java.util.HashSet;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.strtree.STRtree;

public class TrackManager {

	private STRtree index;

	public TrackManager() {
		this.index = new STRtree();
	}

	public void addTrack(Track track) {
		for (TrackPoint tp : track.getPoints())
			index.insert(new Envelope(tp.getCoordinate()), track);
	}

	public ArrayList<Track> queryTracks(Envelope searchWindow) {
		@SuppressWarnings("unchecked")
		HashSet<Track> trac = new HashSet<>(index.query(searchWindow));
		return new ArrayList<>(trac);
	}

}
