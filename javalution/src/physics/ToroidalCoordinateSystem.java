package physics;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import general.Distribution;
import general.Pair;

/**
 * A CoordinateSystem with a toroidal topology.
 * That is, a square grid where each side is connected to its opposite.
 * A cell's neighbors are its Moore neighborhood (up-down-left-right and diagonals).
 */
public class ToroidalCoordinateSystem implements CoordinateSystem {
	
	public ToroidalCoordinateSystem(int gridSize) {
		this.gridSize = gridSize;
		for (int x=0; x<gridSize; x++) {
			this.map.add(new ArrayList<>());
			for (int y=0; y<gridSize; y++) {
				ToroidalCoordinate co = this.new ToroidalCoordinate(x,y);
				this.map.get(x).add(co);
				this.all.add(co); 
			}
		}
		this.randomizer = Distribution.uniform(all);
	}
	
	private final int gridSize;
	private final Set<Coordinate> all = new HashSet<>();
	private final List<List<ToroidalCoordinate>> map = new ArrayList<>();
	private final Distribution<Coordinate> randomizer;

	@Override
	public Coordinate random() {
		return this.randomizer.getRandom();
	}
	
	@Override
	public Coordinate random(Predicate<Coordinate> filter) {
		return Distribution.uniform(CoordinateSystem.filter(this.all, filter)).getRandom();
	}
	
	public ToroidalCoordinate get(int x, int y) {
		return this.map
				.get(((x%this.gridSize)+this.gridSize)%this.gridSize)
				.get(((y%this.gridSize)+this.gridSize)%this.gridSize);
	}

	public Set<Coordinate> neighbors(ToroidalCoordinate middle, int distance) {
		if (middle == null) throw new IllegalArgumentException("null does not have neighbors");
		if (distance < 0) throw new IllegalArgumentException("Cannot find neighbors with negative distance");
		Set<Coordinate> result = new HashSet<>();
		for (int i=-(distance%this.gridSize); i<=distance%this.gridSize; i++)
		for (int j=-distance%this.gridSize; j<=distance%this.gridSize; j++)
			if (i!=0 || j!=0)
				result.add(this.get(middle.left+i,middle.right+j));
		return result;
	}
	
	public class ToroidalCoordinate extends Pair<Integer,Integer> implements CoordinateSystem.Coordinate {
		
		private ToroidalCoordinate(int x, int y) {
			super(x,y);
		}

		@Override
		public String toString() {
			return "("+this.left+","+this.right+")";
		}

		@Override
		public Set<CoordinateSystem.Coordinate> neighbors(int i) {
			return ToroidalCoordinateSystem.this.neighbors(this, i);
		}
	}

	@Override
	public boolean has(Coordinate c) {
		return this.all.contains(c);
	}

}
