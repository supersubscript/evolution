package settings;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import organisms.Organism;
import physics.World;

/**
 * Functional interface representing a fitness function.
 * This is essentially a ToDoubleFunction, but can also rescale the fitness values depending on other factors in the World.
 * For example, normalizing fitness to sum to 1, decreasing fitness if high-fitness neighbors are nearby, etc.
 */
public interface FitnessFunction<O extends Organism<O>> extends Function<O,Double> {
	
	/**
	 * Apply this fitness function directly.
	 * This function does not take information on the Organism's environment and is thus not modified.
	 */
	@Override
	Double apply(O value);
	
	/**
	 * Apply this FitnessFunction to all Organisms in a World.
	 * Unlike apply(O), this function normalizes the fitness values.
	 */
	public default Map<O,Double> apply(World<O> w) {
		return w.getPopulation().stream().collect(Collectors.toMap(Function.identity(), this));
	}
	
	/**
	 * Get the total, non-normalized fitness of the population.
	 */
	public default double totalFitness(World<O> w) {
		return w.getPopulation().stream().mapToDouble(x -> this.apply(x)).sum();
	}
		
	/**
	 * Rescaled FitnessFunction that normalizes all fitnesses in a given World so that the total fitness of all organisms equals 1.
	 */
	public static <O extends Organism<O>> FitnessFunction<O> normalize(FitnessFunction<O> fitness) {
		return new FitnessFunction<O>() {
			@Override
			public Double apply(O value) {
				return fitness.apply(value);
			}
			@Override
			public Map<O, Double> apply(World<O> w) {
				Map<O,Double> map = fitness.apply(w);
				Double totalFitness = map.values().stream().mapToDouble(Double::doubleValue).sum();
				for (Map.Entry<O, Double> entry : map.entrySet()) entry.setValue(entry.getValue() / totalFitness);
				return map;
			}
		};
	}
	
	/**
	 * Rescaled FitnessFunction whose fitness values are a power of the original (typically after normalization).
	 */
	public static <O extends Organism<O>> FitnessFunction<O> power(FitnessFunction<O> fitness, Double p) {
		return new FitnessFunction<O>() {
			@Override
			public Double apply(O value) {
				return fitness.apply(value);
			}
			@Override
			public Map<O, Double> apply(World<O> w) {
				Map<O,Double> map = fitness.apply(w);
				for (Map.Entry<O, Double> entry : map.entrySet()) entry.setValue(Math.pow(entry.getValue(),p));
				return map;
			}
		};
	}
	
	/**
	 * Rescaled FitnessFunction that sets the fitness value of the lowth percentile individual to 0, and the highth percentile individual to 1, linearly interpolating in between.  
	 */
	public static <O extends Organism<O>> FitnessFunction<O> cutoff(FitnessFunction<O> fitness, Double low, Double high, boolean allowHigherThan1) {
		double maxfit = allowHigherThan1 ? Double.MAX_VALUE : 1.;
		return new FitnessFunction<O>() {
			@Override
			public Double apply(O value) {
				return fitness.apply(value);
			}
			@Override
			public Map<O, Double> apply(World<O> w) {
				Map<O,Double> map = fitness.apply(w);
				List<Double> fs = new ArrayList<>(map.values());
				Collections.sort(fs);
				Double lowthP = fs.get((int) Math.round(low * w.getPopulation().size()));
				Double highthP = fs.get((int) Math.round(high * w.getPopulation().size() - 1 ));
				for (Map.Entry<O, Double> entry : map.entrySet()) entry.setValue(Math.max(0, Math.min(maxfit, (entry.getValue() - lowthP)/(highthP - lowthP))));
				return map;
			}
		};
	}

}
