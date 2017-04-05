package settings;

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import general.Distribution;
import organisms.Organism;
import physics.World;

/**
 * Strategy interface for selecting subsets of Organisms from a population, such as when deciding which organism dies or splits.
 */
@FunctionalInterface
public interface SelectionRule<O extends Organism<O>> extends Function<World<O>,Collection<O>> {
	
	/**
	 * Alias of {@link #select(World)}.
	 */
	@Override
	default Collection<O> apply(World<O> t) {
		return this.select(t);
	}

	/**
	 * Which Organisms to select, out of all the Organisms in the World.
	 */
	public Collection<O> select(World<O> w);
	
	/**
	 * Whether one specific organism is selected.
	 */
	public default boolean selectQ(World<O> w, O o) {
		return this.select(w).contains(o);
	}
	
	/**
	 * A selection rule that selects all Organisms.
	 */
	public static <O extends Organism<O>> SelectionRule<O> all() {
		return new SelectionRule<O>() {
			@Override
			public boolean selectQ(World<O> w, O o) {
				return true;
			}
			@Override
			public Set<O> select(World<O> w) {
				return w.getPopulation();
			}
		};
	}
	
	/**
	 * A selection rule that selects no Organisms.
	 */
	public static <O extends Organism<O>> SelectionRule<O> none() {
		return new SelectionRule<O>() {
			@Override
			public boolean selectQ(World<O> w, O o) {
				return false;
			}
			@Override
			public Set<O> select(World<O> w) {
				return new HashSet<O>();
			}
		};
	}
	
	/**
	 * A rule that selects random Organisms.
	 * The probability that any Organism will be selected is exactly rate.
	 */
	public static <O extends Organism<O>> SelectionRule<O> constantRate(double rate) {
		return new SelectionRule<O>() {
			Random rand = new Random();
			@Override
			public boolean selectQ(World<O> w, O o) {
				return rate > this.rand.nextDouble();
			}
			@Override
			public Set<O> select(World<O> w) {
				Set<O> result = new HashSet<>();
				for (O o : w.getPopulation()) {
					if (selectQ(w, o)) result.add(o);
				}
				return result; // TODO make more efficient
			}
		};
	}
	
	/**
	 * A rule that selects random Organisms.
	 * The probability of being selected goes up proportionally with population size, i.e. probability of being selected is
	 * 		rate + populationsize*raterate
	 * Selected Organisms are removed from the pool so that probability of selection never exceeds 1.
	 */
	public static <O extends Organism<O>> SelectionRule<O> proportionalRate(double rate, double raterate) {
		return new SelectionRule<O>() {
			Random rand = new Random();
			@Override
			public Set<O> select(World<O> w) {
				Set<O> result = new HashSet<>();
				for (O o : w.getPopulation()) {
					if (rate + (w.getPopulation().size() - result.size()) * raterate > rand.nextDouble()) result.add(o);
				}
				return result;
			}
		};
	}
	
	/**
	 * Select a random sample from the population.
	 */
	public static <O extends Organism<O>> SelectionRule<O> selectProportion(double proportion) {
		double internalProportion = Math.max(proportion, 0.);
		return new SelectionRule<O>() {
			@Override
			public Collection<O> select(World<O> w) {
				Distribution<O> randomizer = Distribution.uniform(w.getPopulation());
				return randomizer.getRandomWithoutReplacement((int)Math.floor(internalProportion * w.getPopulation().size()));
			}
			@Override
			public String toString() {
				return "Selection rule (proportion "+proportion+")";
			}
		};
	}

	/**
	 * A selection rule that selects a random surplus of organisms, if the population is bigger than the target size.
	 */
	public static <O extends Organism<O>> SelectionRule<O> cullToNumber(int popSize) {
		return new SelectionRule<O>() {
			@Override
			public Collection<O> select(World<O> w) {
				Distribution<O> randomizer = Distribution.uniform(w.getPopulation());
				return randomizer.getRandom(Math.max(0, w.getPopulation().size() - popSize));
			}
			@Override
			public String toString() {
				return "Selection rule (culler "+popSize+")";
			}
		};
	}

	/**
	 * Represents a mechanism based on the number of neighbors.
	 * All Organisms with i or less neighbors are selected.
	 */
	public static <O extends Organism<O>> SelectionRule<O> maxNeighbors(int i, int dist) {
		return new SelectionRule<O>() {
			@Override
			public Set<O> select(World<O> w) {
				LocationFilter fullTest = LocationFilter.full();
				Set<O> result = new HashSet<>();
				for (O o : w.getPopulation()) {
					int neighbors = fullTest.count(w, w.getLocation(o).neighbors(dist));
					if (neighbors <= i) result.add(o); 
				}
				return result;
			}
		};
	}

	/**
	 * Represents a mechanism based on the number of neighbors.
	 * All Organisms with i or more neighbors are selected.
	 */
	public static <O extends Organism<O>> SelectionRule<O> minNeighbors(int i, int dist) {
		return new SelectionRule<O>() {
			@Override
			public Set<O> select(World<O> w) {
				LocationFilter fullTest = LocationFilter.full();
				Set<O> result = new HashSet<>();
				for (O o : w.getPopulation()) {
					int neighbors = fullTest.count(w, w.getLocation(o).neighbors(dist));
					if (neighbors >= i) result.add(o); 
				}
				return result;
			}
		};
	}

	/**
	 * Selection rule where the probability of selecting any individual x in the population is equal to a*f(x)^p + b.
	 */
	public static <O extends Organism<O>> SelectionRule<O> rouletteWheelSelection(FitnessFunction<O> fitness) {
		return new SelectionRule<O>() {
			Random rand = new Random();
			@Override
			public Set<O> select(World<O> w) {
				Map<O,Double> fs = fitness.apply(w);
				return fs.entrySet().stream()
						.filter(x -> x.getValue() > rand.nextDouble())
						.map(Map.Entry::getKey)
						.collect(Collectors.toSet());
			}
		};
	}
	
	/**
	 * Select a random sample from the population.
	 * Probability of being selected is proportional to fitness function.
	 * Random selection is done without replacement.
	 */
	public static <O extends Organism<O>> SelectionRule<O> selectProportion(FitnessFunction<O> fitness, double proportion) {
		double internalProportion = Math.max(proportion, 0.);
		return new SelectionRule<O>() {
			@Override
			public Collection<O> select(World<O> w) {
				Distribution<O> randomizer = Distribution.weighted(fitness.apply(w));
				return randomizer.getRandomWithoutReplacement((int)Math.floor(internalProportion * w.getPopulation().size()));
			}
			@Override
			public String toString() {
				return "Selection rule (proportion "+proportion+"; fitness "+fitness+")";
			}
		};
	}
	
	/**
	 * A selection rule that selects a random surplus of organisms, if the population is bigger than the target size.
	 * Probability of being selected is proportional to fitness.
	 */
	public static <O extends Organism<O>> SelectionRule<O> cullToNumber(FitnessFunction<O> fitness, int popSize) {
		return new SelectionRule<O>() {
			@Override
			public Collection<O> select(World<O> w) {
				Distribution<O> randomizer = Distribution.weighted(fitness.apply(w));
				return randomizer.getRandom(Math.max(0, w.getPopulation().size() - popSize));
			}
			@Override
			public String toString() {
				return "Selection rule (culler "+popSize+")";
			}
		};
	}
}