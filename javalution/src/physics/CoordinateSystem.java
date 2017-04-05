package physics;

import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import general.Distribution;

/**
 * This class serves to pinpoint a specific way of dealing with Coordinates.
 * Implementing classes should also provide an implementation of Coordinates.
 * If the size of the world is finite, implementations can also serve as a flyweight for those Coordinates.
 *
 */
public interface CoordinateSystem {
	
	/**
	 * Filter a set of Coordinates.
	 */
	public static <C extends Coordinate> Set<C> filter(Set<C> coordinates, Predicate<? super C> filter) {
		return coordinates.stream().filter(filter).collect(Collectors.toSet());
	}
	
	/**
	 * Check if a specified Coordinate belongs to this CoordinateSystem.
	 */
	public abstract boolean has(Coordinate c);
	
	/**
	 * Get a random Coordinate from this System.
	 */
	public abstract Coordinate random();
	
	/**
	 * Get a random Coordinate from this System that does not have 
	 */
	public abstract Coordinate random(Predicate<Coordinate> filter);
	
	/**
	 * Interface representing a Coordinate in a certain CoordinateSystem.
	 * The Coordinate interface provides only neighborhood structure.
	 * Any other properties, such as number of possible Coordinates or symmetry are encapsulated.
	 */
	public interface Coordinate {
		public Set<Coordinate> neighbors(int i);
		public default Set<Coordinate> neighbors() {
			return this.neighbors(1);
		}
		public default Coordinate randomNeighbor(int i) {
			return Distribution.uniform(this.neighbors(i)).getRandom();
		}
		public default Coordinate randomNeighbor() {
			return this.randomNeighbor(1);
		}
		public default Coordinate randomNeighbor(int i, Predicate<Coordinate> filter) {
			return Distribution.uniform(CoordinateSystem.filter(this.neighbors(i), filter)).getRandom();
		}
		public default Coordinate randomNeighbor(Predicate<Coordinate> filter) {
			return this.randomNeighbor(1, filter);
		}
	}

}
