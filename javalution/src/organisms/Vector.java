package organisms;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import general.Distribution;
import general.Sequence;
import settings.CrossOverOperator;
import settings.MutationOperator;
import settings.RealReader;

/**
 * Organism that represents a variable amount of real numbers.
 * The real numbers are defined in the same way as for the Number class,
 * except that the genome is first split into segments separated by a given Sequence.
 * Each segment is then interpreted as a separate real number. 
 * If no separation into multiple numbers is required, set separator to null. Then this class behaves like Number.
 */
public class Vector extends Organism<Vector> {
	
	private Vector(Sequence genome, List<Double> phi) {
		this.genome = genome;
		this.phi = phi;
	}
	
	private final Sequence genome;
	private final List<Double> phi;
	
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
		builder.append("          <---------->          ");
		builder.append(this.genome.toString());
		builder.append("   ////////   ");
		builder.append(this.genome.split(Sequence.create(new boolean[]{false, true, true, true, true, true, false})));
		return builder.toString();
	}
	
	public static class Factory implements Organism.Factory<Vector> {
		
		/**
		 * Create a new Vector Factory.
		 * @param mut
		 * 		Mutation operator that is called whenever an Organism multiplies asexually or sexually.
		 * @param initialGenomeSize
		 * 		Random new Vectors (Factory.random()) will have random genomes with a length gotten from this distribution.
		 * @param reader
		 * 		RealReader that describes how to read real numbers from this Genome.
		 * @param crossover
		 * 		CrossOverOperator that describes how to combine two genomes into a new Vector.
		 * @param separator
		 * 		Which sequence separates the constituents of a Vector.
		 * 		If no separation is wanted, use separator=null. Then this class behaves similarly to a Number.
		 */
		public Factory(
				MutationOperator mut,
				Distribution<Integer> initialGenomeSize,
				RealReader reader,
				CrossOverOperator crossover,
				Sequence separator) {
			
			if (mut == null) throw new IllegalArgumentException("null cannot be a mutation operator. Use MutationOperator.trivial() instead.");
			if (initialGenomeSize == null) throw new IllegalArgumentException("null is not a valid initial genome size distribution");
			if (reader == null) throw new IllegalArgumentException("null is not a valid RealReader");
			if (crossover == null) throw new IllegalArgumentException("null is nto a valid CrossOverOperator");
			
			this.mut = mut;
			this.initialGenomeSize = initialGenomeSize;
			this.reader = reader;
			this.crossover = crossover;
			this.separator = separator;
		}
		
		private final MutationOperator mut;
		private final Distribution<Integer> initialGenomeSize;
		private final RealReader reader;
		private final CrossOverOperator crossover;
		private final Sequence separator;
		
		private List<Double> read(Sequence seq) {
			List<Double> result = new ArrayList<>();
			for (Sequence s : seq.split(this.separator)) {
				result.add(this.reader.apply(s));
			}
			return Collections.unmodifiableList(result);
		}

		@Override
		public Vector random() {
			Sequence     seq = Sequence.random(this.initialGenomeSize.get());
			List<Double> phi = this.read(seq);
			return new Vector(seq, phi);
		}

		@Override
		public Vector copy(Vector o) {
			Sequence     seq = o.genome;
			List<Double> phi = this.read(seq);
			return new Vector(seq, phi);
		}

		@Override
		public Vector split(Vector o) {
			Sequence     seq = this.mut.mutate(o.genome);
			List<Double> phi = this.read(seq);
			return new Vector(seq, phi);
		}
		
		@Override
		public Vector sex(Vector mommy, Vector daddy) {
			Sequence     seq = this.mut.mutate(this.crossover.apply(mommy.genome, daddy.genome));
			List<Double> phi = this.read(seq);
			return new Vector(seq, phi);
		}
		
	}

}
