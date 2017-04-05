package organisms;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.function.ToDoubleFunction;

import general.Distribution;
import general.Sequence;
import physics.World;
import settings.CrossOverOperator;
import settings.MutationOperator;
import settings.RealReader;
import settings.SelectionRule;

/**
 * A class much like Vector but with an inner state.
 * A Plant takes up light that is stored as energy.
 * Energy is depleted when reproducing (cost is shared between parents for sexual reproduction).
 * Reproduction only continues when energy is above zero.
 * Plant.tresholdEnergy() also implements a good DeathRule for use with this class; it kills all Plants with energy below a given treshold.
 */
public class Plant extends Organism<Plant> {
	
	private Plant(Sequence genome, List<Double> phi, double energy) {
		this.genome = genome;
		this.phi = phi;
		this.energy = energy;
	}
	
	private final Sequence genome;
	private final List<Double> phi;
	private double energy;
	
	public List<Double> getValues() {
		return this.phi;
	}
	
	public double getValue(int index, double def) {
		if (this.phi.size() <= index) return def;
		return this.phi.get(index);
	}
	
	public double getValue(int index) {
		return this.getValue(index, 0.);
	}
	
	public double getEnergy() {
		return this.energy;
	}
	
	@Override
	public void tick(World<Plant> w) {
		super.tick(w);
		this.energy += w.getLight(this);
	}
	
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		DecimalFormat format = new DecimalFormat("###.##");
		String sep = "[";
		for (double d : this.getValues()) {
			builder.append(sep+format.format(d));
			sep = ", ";
		}
		builder.append("]");
		builder.append(" (" + format.format(this.getEnergy()) + ")");
		builder.append("          <---------->          ");
		builder.append(this.genome.toString());
		return builder.toString();
	}
	
	public static class Factory implements Organism.Factory<Plant> {
		
		public Factory(
				Distribution<Integer> initialGenomeSize,
				MutationOperator mut,
				CrossOverOperator crossover,
				RealReader realReader,
				Sequence separator,
				double initialEnergy,
				ToDoubleFunction<Plant> reproductionCost
				) {
			
			this.initialGenomeSize = initialGenomeSize;
			this.mut = mut;
			this.crossover = crossover;
			this.realReader = realReader;
			this.separator = separator;
			this.initialEnergy = initialEnergy;
			this.reproductionCost = reproductionCost;
			
		}
		
		private final Distribution<Integer> initialGenomeSize;
		private final MutationOperator mut;
		private final CrossOverOperator crossover;
		private final RealReader realReader;
		private final Sequence separator;
		private final double initialEnergy;
		private final ToDoubleFunction<Plant> reproductionCost;
		
		private List<Double> read(Sequence seq) {
			List<Double> result = new ArrayList<>();
			for (Sequence s : seq.split(this.separator)) {
				result.add(this.realReader.apply(s));
			}
			return Collections.unmodifiableList(result);
		}


		@Override
		public Plant random() {
			Sequence     seq = Sequence.random(this.initialGenomeSize.get());
			List<Double> phi = this.read(seq);
			return new Plant(seq, phi, initialEnergy);
		}

		@Override
		public Plant copy(Plant o) {
			return new Plant(o.genome, o.phi, o.energy);
		}

		@Override
		public Plant split(Plant o) {
			o.energy        -= this.reproductionCost.applyAsDouble(o);
			if (o.energy <= 0.) return null;
			Sequence seq     = this.mut.mutate(o.genome);
			List<Double> phi = this.read(seq);
			return new Plant(seq, phi, initialEnergy);
		}

		@Override
		public Plant sex(Plant mommy, Plant daddy) {
			mommy.energy    -= this.reproductionCost.applyAsDouble(mommy)/2.;
			daddy.energy    -= this.reproductionCost.applyAsDouble(daddy)/2.;
			if (mommy.energy <= 0. || daddy.energy <= 0.) return null;
			Sequence seq     = this.mut.mutate(this.crossover.apply(mommy.genome, daddy.genome));
			List<Double> phi = this.read(seq);
			return new Plant(seq, phi, initialEnergy);
		}
		
	}
	
	/**
	 * A SelectionRule that selects Plants with energy below a certain treshold value.
	 */
	public static SelectionRule<Plant> tresholdEnergy(double treshold) {
		return new SelectionRule<Plant>() {
			@Override
			public Collection<Plant> select(World<Plant> w) {
				Collection<Plant> result = new HashSet<>();
				for (Plant p : w.getPopulation())
					if (p.getEnergy() < treshold)
						result.add(p);
				return result;
			}
		};
	}

}
