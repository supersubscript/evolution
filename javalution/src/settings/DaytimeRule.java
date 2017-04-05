package settings;

import general.Distribution;
import organisms.Organism;
import physics.World;

/**
 * A functional interface representing the functionality of calculating the amount of light present at a particular tick number.
 */
@FunctionalInterface
public interface DaytimeRule<O extends Organism<O>> {
	
	public double light(World<O> w, O o);
	
	/**
	 * A trivial implementation with no light.
	 */
	public static <O extends Organism<O>> DaytimeRule<O> trivial() {
		return DaytimeRule.constant(0.);
	}
	
	/**
	 * The light that falls is constant.
	 */
	public static <O extends Organism<O>> DaytimeRule<O> constant(double constant) {
		return new DaytimeRule<O>() {
			@Override
			public double light(World<O> w, O o) {
				return constant;
			}
		};
	}
	
	/**
	 * Create a DaytimeRule rule where light is distributed evenly along all Organisms.
	 */
	public static <O extends Organism<O>> DaytimeRule<O> distribute(DaytimeRule<O> rule) {
		return new DaytimeRule<O>() {
			@Override
			public double light(World<O> w, O o) {
				return rule.light(w, o) / w.getPopulation().size();
			}
		};
	}
	
	/**
	 * Create a DaytimeRule rule with noise.
	 */
	public static <O extends Organism<O>> DaytimeRule<O> noise(DaytimeRule<O> rule, Distribution<Double> noise) {
		return new DaytimeRule<O>() {
			@Override
			public double light(World<O> w, O o) {
				return rule.light(w, o) + noise.getRandom();
			}
		};
	}

}
