package settings;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.function.Predicate;

import physics.CoordinateSystem;
import physics.CoordinateSystem.Coordinate;
import physics.Habitat;

/**
 * Represents a filter on locations.
 */
@FunctionalInterface
public interface LocationFilter {

	public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate);
	
	public default Set<CoordinateSystem.Coordinate> filter(Habitat<?> w, Set<CoordinateSystem.Coordinate> coordinates) {
		Set<CoordinateSystem.Coordinate> result = new HashSet<>();
		for (CoordinateSystem.Coordinate co : coordinates) {
			if (this.test(w, co)) result.add(co);
		}
		return result;
	}
	
	public default List<CoordinateSystem.Coordinate> filter(Habitat<?> w, List<CoordinateSystem.Coordinate> coordinates) {
		List<CoordinateSystem.Coordinate> result = new ArrayList<>();
		for (CoordinateSystem.Coordinate co : coordinates) {
			if (this.test(w, co)) result.add(co);
		}
		return result;
	}
	
	public default int count(Habitat<?> w, Collection<CoordinateSystem.Coordinate> coordinates) {
		int result = 0;
		for (CoordinateSystem.Coordinate co : coordinates) if (this.test(w, co)) result++;
		return result;
	}
	
	public default Predicate<Coordinate> concretize(Habitat<?> w) {
		return (c) -> this.test(w, c);
	}
	
	/**
	 * A filter that is trivially true for any argument.
	 */
	public static LocationFilter trivialTrue() {
		return new LocationFilter() {
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				return true;
			}
			@Override
			public int count(Habitat<?> w, Collection<CoordinateSystem.Coordinate> coordinates) {
				return coordinates.size();
			}
			@Override
			public List<CoordinateSystem.Coordinate> filter(Habitat<?> w, List<CoordinateSystem.Coordinate> coordinates) {
				return coordinates;
			}
			@Override
			public Set<CoordinateSystem.Coordinate> filter(Habitat<?> w, Set<CoordinateSystem.Coordinate> coordinates) {
				return coordinates;
			}
		};
	}
	/**
	 * A filter that is trivially false for any argument.
	 */
	public static LocationFilter trivialFalse() {
		return new LocationFilter() {
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				return false;
			}
			@Override
			public int count(Habitat<?> w, Collection<CoordinateSystem.Coordinate> coordinates) {
				return 0;
			}
			@Override
			public List<CoordinateSystem.Coordinate> filter(Habitat<?> w, List<CoordinateSystem.Coordinate> coordinates) {
				return new ArrayList<>();
			}
			@Override
			public Set<CoordinateSystem.Coordinate> filter(Habitat<?> w, Set<CoordinateSystem.Coordinate> coordinates) {
				return new HashSet<>();
			}
		};
	}
	/**
	 * A filter that is randomly true or false with a certain probability, regardless of the arguments.
	 */
	public static LocationFilter random(double prob) {
		return new LocationFilter() {
			Random rand = new Random();
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				return prob > this.rand.nextDouble();
			}
		};
	}
	
	/**
	 * A filter that is true iff the location contains an Organism.
	 */
	public static LocationFilter full() {
		return new LocationFilter() {
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				return w.isOccupied(coordinate);
			}
		};
	}
	
	/**
	 * A filter that is true iff the location does not contain an Organism.
	 */
	public static LocationFilter empty() {
		return new LocationFilter() {
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				return !w.isOccupied(coordinate);
			}
		};
	}

	/**
	 * Represents a filter based on the number of neighbors that are filled with Organisms.
	 * All locations with i or less neighbors are selected.
	 */
	public static LocationFilter maxNeighbors(int i, int dist) {
		return new LocationFilter() {
			LocationFilter testOccupied = LocationFilter.full();
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				int filled = testOccupied.count(w, coordinate.neighbors(dist));
				return filled <= i;
			}
		};
	}
	
	/**
	 * Represents a filter based on the number of neighbors that are filled with Organisms.
	 * All locations with i or more neighbors are selected.
	 */
	public static LocationFilter minNeighbors(int i, int dist) {
		return new LocationFilter() {
			LocationFilter testOccupied = LocationFilter.full();
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				int filled = testOccupied.count(w, coordinate.neighbors(dist));
				return filled >= i;
			}
		};
	}
	
	/**
	 * Represents a filter based on the number of neighbors that are filled with Organisms.
	 * All locations with exactly i neighbors is selected.
	 */
	public static LocationFilter neighbors(int i, int dist) {
		return new LocationFilter() {
			LocationFilter testOccupied = LocationFilter.full();
			@Override
			public boolean test(Habitat<?> w, CoordinateSystem.Coordinate coordinate) {
				int filled = testOccupied.count(w, coordinate.neighbors(dist));
				return filled == i;
			}
		};
	}
}