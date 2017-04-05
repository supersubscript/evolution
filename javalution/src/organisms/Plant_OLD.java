package organisms;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.function.Function;

import general.Distribution;
import general.Sequence;
import physics.World;
import settings.MutationOperator;
import settings.SelectionRule;

/**
 * An implementation of the Plant_OLD proof-of-concept.
 * A Plant_OLD's genome consists of a list of doubles, each representing the probability of replicating at a particular moment in time.
 * The Plant_OLD's internal state consists of an energy level. At each tick, the World's photo-energy is divided among all Plants.
 * New plants are generated with the minimal treshold (thus, in times of negative net energy income, reproduction always fails)
 * Design requirement: Plants must be put in a SeasonWorld in order to receive energy.
 * Design requirement: The World Settings must contain a rule 
 */
public class Plant_OLD extends Organism<Plant_OLD> {
	
	private Plant_OLD(Sequence genome, List<Double> phenotype, double energyTreshold) {
		this.genome = genome;
		this.phenotype = phenotype;
		this.energyTreshold = energyTreshold;
	}
	
	private final Sequence genome;
	private final List<Double> phenotype;
	private double energy;
	private final double energyTreshold;
	
	public double divisionProbability(double time) {
		if (this.phenotype.size() == 0) return 0.;
		int index = (int) Math.floor(time * (this.phenotype.size() - Double.MIN_VALUE));
		return this.phenotype.get(index);
	}
	
	@Override
	public void tick(World<Plant_OLD> w) {
		// Absorb energy
		this.energy += w.getLight(this) / (w.getPopulation().size());
	}

	public static class Factory implements Organism.Factory<Plant_OLD> {
		public Factory(
				MutationOperator mut,
				Distribution<Integer> initialGenomeSize,
				Sequence geneStart,
				Sequence geneStop,
				Function<Sequence,Double> realReader,
				double energyTreshold,
				double initialEnergy
				) {
			this.mut = mut;
			this.initialGenomeSize = initialGenomeSize;
			this.geneStart = geneStart;
			this.geneStop = geneStop;
			this.realReader = realReader;
			this.energyTreshold = energyTreshold;
			this.initialEnergy = initialEnergy;
		}
		
		private final MutationOperator mut;
		private final Distribution<Integer> initialGenomeSize;
		private final Sequence geneStart;
		private final Sequence geneStop;
		private final Function<Sequence,Double> realReader;
		private final double initialEnergy;
		private final double energyTreshold;

		@Override
		public Plant_OLD random() {
			Sequence genome = Sequence.random(initialGenomeSize.getRandom());
			List<Double> phenotype = new ArrayList<Double>();
			for (Sequence s : genome.genes(geneStart, geneStop)) phenotype.add(realReader.apply(s));
			Plant_OLD result = new Plant_OLD(genome, phenotype, energyTreshold);
			result.energy = initialEnergy;
			return result;
		}

		@Override
		public Plant_OLD copy(Plant_OLD o) {
			Plant_OLD po = o;
			return new Plant_OLD(po.genome, po.phenotype, energyTreshold);
		}

		@Override
		public Plant_OLD split(Plant_OLD o) {
			o.energy = o.energy - energyTreshold;
			if (o.energy < energyTreshold) return null;
			Sequence genome = mut.mutate(o.genome);
			List<Double> phenotype = new ArrayList<Double>();
			for (Sequence s : genome.genes(geneStart, geneStop)) phenotype.add(realReader.apply(s));
			Plant_OLD result = new Plant_OLD(genome, phenotype, energyTreshold);
			o.energy = o.energy / 2;
			result.energy = o.energy;
			return result;
		}

		@Override
		public Plant_OLD sex(Plant_OLD mommy, Plant_OLD daddy) {
			// TODO Auto-generated method stub
			return null;
		}
	}
	
	public String genotypeString() {
		return this.genome.toString();
	}
	
	public String phenotypeString(int digitPrecision) {
		DecimalFormat format = new DecimalFormat("###.###");
		StringBuilder builder = new StringBuilder();
		builder.append('(');
		String prefix = "";
		for (double d : this.phenotype) {
			builder.append(prefix); prefix = "-";
			builder.append(format.format(Math.min(1.,d)));
		}
		builder.append(')');
		return builder.toString();
	}
	
	@Override
	public String toString() {
		return "[" + this.genotypeString() + ", " + this.phenotypeString(1000) + ", " + this.energy + "]";
	}
	
	public static SelectionRule<Plant_OLD> deathRule() {
		return new SelectionRule<Plant_OLD>() {
			@Override
			public Collection<Plant_OLD> select(World<Plant_OLD> w) {
				Collection<Plant_OLD> result = new HashSet<Plant_OLD>();
				for (Plant_OLD p : w.getPopulation())
					if (p.energy < p.energyTreshold)
						result.add(p);
				return result;
			}
		};
	}
	
	public static SelectionRule<Plant_OLD> splitRule() {
		return new SelectionRule<Plant_OLD>() {
			@Override
			public Collection<Plant_OLD> select(World<Plant_OLD> w) {
				Collection<Plant_OLD> result= new HashSet<Plant_OLD>();
				for (Plant_OLD p : w.getPopulation())
					if (p.energy > p.energyTreshold && w.settings.rand.nextDouble() < p.divisionProbability(w.getTime()))
						result.add(p);
				return result;
			}
		};
	}
	
}
