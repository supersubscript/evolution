package settings;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.function.Supplier;

import general.Sequence;

/**
 * Strategy interface for easy access to Functions that read Sequences into real values.
 */
@FunctionalInterface
public interface RealReader extends Function<Sequence,Double> {
	
	/**
	 * Reads a binary sequence in regular base-2 encoding.
	 * Note that the direction of reading is reversed!
	 */
	public static RealReader binary() {
		return new RealReader() {
			@Override
			public Double apply(Sequence t) {
				double result = 0.;
				for (int i=0; i<t.length(); i++)
					if (t.get(i)) result += Math.pow(2, i);
				return result;
			}
		};
	}
	
	/**
	 * Reads a binary sequence using the regular "binary reflected" Gray encoding.
	 * Note that the direction of reading is reversed!
	 */
	public static RealReader gray() {
		return new RealReader() {
			RealReader binaryReader = RealReader.binary();
			@Override
			public Double apply(Sequence t) {
				int bin = binaryReader.apply(t).intValue();
				int result = bin;
				while ((bin >>>= 1) != 0) // Unintelligable code pulled from internet.
					result ^= bin;
				return (double)result;
			}
		};
	}

	/**
	 * This RealReader reads a sequence in a multiplicative manner.
	 * Each position in the sequence is assigned a random (but static) weight, calculated according to the given weightDistribution. (Alternatively a different Supplier can be used to get deterministic weights)
	 * These weights are summed over the sequence, where each term gets a positive sign if its corresponding bit is 1, or negative if the bit is 0.
	 * The result is the exponential of this sum.
	 * Wider weightDistributions lead to more range and less accuracy/stability.
	 */
	public static RealReader multiplicative(Supplier<Double> weightDistribution) {
		return new RealReader() {
			private final List<Double> weights = new ArrayList<>();
			
			private double weight(int i) {
				while (weights.size() <= i) {
					double newWeight = weightDistribution.get();
					weights.add(newWeight);
				}
				return weights.get(i);
			}
			
			@Override
			public Double apply(Sequence t) {
				double result = 0;
				for (int i=0; i<t.length(); i++) {
					result += t.get(i) ? weight(i) : -weight(i);
				}
//				Function<Double,Double> fit = (d) -> {
//					double val = Math.max(d, 0.);
//					return val > 2.
//							? 0.
//							: val * Math.sin(10. * Math.PI * val) + 2.;
//				};
//				DecimalFormat format = new DecimalFormat("###.###");
//				if (Math.random() < 0.01); System.out.println(""
//						+ t + " -- "
//						+ format.format(result) + " -- "
//						+ format.format(Math.exp(result)) + " -- "
//						+ format.format(fit.apply(Math.exp(result)))
//						);
				return Math.exp(result);
			}
		};
	}
	
	/**
	 * A RealReader strategy that recognizes certain subsequences in a sequence.
	 * For each such subsequence, the result is multiplied by a given value. The value depends on the region surrounding the recognition sequence.
	 * In addition, inversions of the consensus sequence are also recognized, and result is divided by the corresponding value.
	 * @param consensus
	 * 			The consensus region signaling numeric factors.
	 * @param valueDiameter
	 * 			The amount of bits that should be taken to the left and right of each consensus sequence to calculate the corresponding factor value.
	 * @param reader
	 * 			The reader that converts the bits surrounding each consensus sequence to a factor value.
	 */
	public static RealReader consensus(Sequence consensus, int valueDiameter, RealReader reader) {
		if (valueDiameter < 0) throw new IllegalArgumentException("Cannot use negative value diameter");
		if (reader == null) throw new IllegalArgumentException("Cannot use null reader");
		if (consensus == null || consensus.length() == 0) throw new IllegalArgumentException("Cannot use null or empty consensus sequence");
		return new RealReader() {
			
			@Override
			public Double apply(Sequence t) {
				double result = 1.;
				for (int start : t.findAll(consensus, true)) {
					if (start < valueDiameter || start > t.length() - consensus.length() - valueDiameter) continue;
					Sequence surrounding = Sequence.cat(
							t.get(start-valueDiameter, start),
							t.get(start+consensus.length(), start+consensus.length()+valueDiameter)
							); 
//					System.out.println(surrounding);
					result *= reader.apply(surrounding);
				}
				for (int start : t.findAll(consensus.complement(), true)) {
					if (start < valueDiameter || start > t.length() - consensus.length() - valueDiameter) continue;
					Sequence surrounding = Sequence.cat(
							t.get(start-valueDiameter, start),
							t.get(start+consensus.length(), start+consensus.length()+valueDiameter)
							); 
					assert reader.apply(surrounding) != 0.;
//					System.err.println(surrounding);
					result /= reader.apply(surrounding);
				}
				return result;
			}
		};
	}
	
	/**
	 * A RealReader strategy that applies a function to an existing RealReader.
	 * That is, the resulting RealReader returns the same value as the supplied reader, but after passing it to the supplied function.
	 */
	public static RealReader applyFunction(RealReader reader, Function<Double,Double> fun) {
		return new RealReader() {
			@Override
			public Double apply(Sequence t) {
				return fun.apply(reader.apply(t));
			}
		};
	}
	
}
