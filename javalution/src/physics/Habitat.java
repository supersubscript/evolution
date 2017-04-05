package physics;

import java.util.List;
import java.util.Set;

import general.Bijection;
import general.Distribution;
import organisms.Organism;

/**
 * A place where Organisms live; connects Organisms to specific Coordinates in a CoordinateSystem.
 * This class serves to encapsulate location-related methods of World.
 */
public class Habitat<O extends Organism<O>> {
	
	public Habitat(CoordinateSystem cs) {
		this.coordinateSystem = cs;
	}
	
	private final CoordinateSystem coordinateSystem;
	final private Bijection<CoordinateSystem.Coordinate,O> grid = new Bijection<>();
	
	////// Inspectors
	
	public Set<O> getPopulation() {
		return this.grid.bvalues();
	}
	
	public Set<CoordinateSystem.Coordinate> getOccupiedLocations() {
		return this.grid.avalues();
	}
	
	public CoordinateSystem.Coordinate getRandomLocation() {
		return this.coordinateSystem.random();
	}
	
	public CoordinateSystem.Coordinate getRandomEmptyLocation() {
		return this.coordinateSystem.random((co) -> !this.isOccupied(co));
	}
	
	public O getRandomOrganism() {
		return Distribution.uniform(this.getPopulation()).getRandom();
	}
	
	public boolean contains(O o) {
		return this.grid.hasB(o);
	}
	
	public CoordinateSystem.Coordinate getLocation(O o) {
		return this.grid.b2a(o);
	}
	
	public List<CoordinateSystem.Coordinate> getLocation(List<O> os) {
		return this.grid.b2a(os);
	}
	
	public Set<CoordinateSystem.Coordinate> getLocation(Set<O> os) {
		return this.grid.b2a(os);
	}
	
	public boolean isOccupied(CoordinateSystem.Coordinate coordinates) {
		return this.getOccupant(coordinates)!=null;
	}
	
	public O getOccupant(CoordinateSystem.Coordinate coordinate) {
		return this.grid.a2b(coordinate);
	}
	
	public List<O> getOccupant(List<CoordinateSystem.Coordinate> coordinates) {
		return this.grid.a2b(coordinates);
	}
	
	public Set<O> getOccupant(Set<CoordinateSystem.Coordinate> coordinates) {
		return this.grid.a2b(coordinates);
	}
	
	public Set<O> getNeighbors(O o, int distance) {
		return this.getOccupant(this.getLocation(o).neighbors(distance));
	}
	
	public Set<O> getNeighbors(O o) {
		return this.getOccupant(this.getLocation(o).neighbors());
	}
	
	////// Mutators
	
	public void kill(O o) {
		if (o != null) {
			o.die();
			this.grid.removeB(o);
		}
	}
	
	public void kill(CoordinateSystem.Coordinate coordinates) {
		this.kill(this.getOccupant(coordinates));
	}
	
	/**
	 * Spawn a new Organism in this Habitat.
	 */
	public void spawnOrganism(O o, CoordinateSystem.Coordinate coordinate) {
		if (o == null) return;
		if (coordinate == null) return;
		if (this.contains(o)) throw new IllegalArgumentException("Cannot spawn Organism that already exists");
		if (!this.coordinateSystem.has(coordinate)) throw new IllegalArgumentException("Coordinate does not belong to the right CoordinateSystem");
		this.grid.add(coordinate, o);
		o.birth();
	}

}
