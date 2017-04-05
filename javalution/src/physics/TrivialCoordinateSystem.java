package physics;

import java.util.HashSet;
import java.util.Set;
import java.util.function.Predicate;

import general.Distribution;

/**
 * A CoordinateSystem set in a trivial metric space of finite size.
 * That is, every cell connects to every other cell.
 */
public class TrivialCoordinateSystem implements CoordinateSystem {
	
	public TrivialCoordinateSystem(int size) {
		this.all = new HashSet<>();
		for (int i=0; i<size; i++) this.all.add(this.new TrivialCoordinate());
		this.randomizer = Distribution.uniform(this.all);
	};
	
	private final Set<Coordinate> all;
	private final Distribution<Coordinate> randomizer;
	
	@Override
	public boolean has(Coordinate c) {
		return this.all.contains(c);
	}
	
	@Override
	public Coordinate random() {
		return this.randomizer.getRandom();
	}
	
	@Override
	public Coordinate random(Predicate<Coordinate> filter) {
		return Distribution.uniform(CoordinateSystem.filter(this.all, filter)).getRandom();
	}
	
	public class TrivialCoordinate implements CoordinateSystem.Coordinate {
		public TrivialCoordinate() {};
		
		@Override
		public Set<Coordinate> neighbors(int i) {
			if (i <= 0) return new HashSet<>();
			Set<Coordinate> result = new HashSet<>(TrivialCoordinateSystem.this.all);
			result.remove(this);
			return result;
		}
		
	}

}
