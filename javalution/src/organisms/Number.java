package organisms;

import java.text.DecimalFormat;

import general.Distribution;
import general.Sequence;
import settings.CrossOverOperator;
import settings.MutationOperator;
import settings.RealReader;

/**
 * Organism that simply represents a single real number.
 * Experiments can select numbers for fitness functions based on their number value, as returned by Number.value();
 * The way that the Number's genome encodes for the number can be set in the factory settings, as well as the mutation operators.
 */
public class Number extends Organism<Number> {
	
	private Number(Sequence genotype, double phenotype) {
		this.genotype  = genotype;
		this.phenotype = phenotype;
	}
	
	private final Sequence genotype;
	private final double   phenotype;
	
	/**
	 * The value that this Number represents.
	 */
	public double getValue() {
		return this.phenotype;
	}
	
	public static class Factory implements Organism.Factory<Number> {
		
		
		/**
		 * Create a new Number Factory.
		 * @param mut
		 * 		Mutation operator that is called whenever an Organism multiplies asexually or sexually.
		 * @param initialGenomeSize
		 * 		Random new Vectors (Factory.random()) will have random genomes with a length gotten from this distribution.
		 * @param reader
		 * 		RealReader that describes how to read real numbers from this Genome.
		 * @param crossover
		 * 		CrossOverOperator that describes how to combine two genomes into a new Organism.
		 */
		public Factory(
				MutationOperator mut,
				Distribution<Integer> initialGenomeSize,
				RealReader reader,
				CrossOverOperator crossover) {
			
			this.mut = mut;
			this.initialGenomeSize = initialGenomeSize;
			this.reader = reader;
			this.crossover = crossover;
		}
		
		private final MutationOperator mut;
		private final Distribution<Integer> initialGenomeSize;
		private final RealReader reader;
		private final CrossOverOperator crossover;
		
		@Override
		public Number random() {
			Sequence seq = Sequence.random(this.initialGenomeSize.getRandom());
			double   phi = this.reader.apply(seq);
			return new Number(seq, phi);
		}
		@Override
		public Number copy(Number o) {
			return new Number(o.genotype.clone(), o.phenotype);
		}
		@Override
		public Number split(Number o) {
			Sequence seq = mut.mutate(o.genotype);
			double   phi = this.reader.apply(seq);
			return new Number(seq, phi);
		}
		@Override
		public Number sex(Number mommy, Number daddy) {
			Sequence seq = this.crossover.apply(mommy.genotype, daddy.genotype);
			double   phi = this.reader.apply(seq);
			return new Number(seq, phi);
		}	
	}
	
	@Override
	public String toString() {
		DecimalFormat format = new DecimalFormat("###.##");
		return this.genotype + " --> " + format.format(this.getValue());
	}

}
