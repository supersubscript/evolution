package settings;

import java.util.Random;
import java.util.function.BiFunction;

import general.Pair;
import general.Sequence;

/**
 * Strategy interface that represents the crossover operator of a digital population.
 */
@FunctionalInterface
public interface CrossOverOperator extends BiFunction<Sequence,Sequence,Sequence> {
	
	/**
	 * CrossOverOperator that simply returns a random parent without crossover.
	 */
	public static CrossOverOperator none() {
		return new CrossOverOperator() {
			Random rand = new Random();
			@Override
			public Sequence apply(Sequence t, Sequence u) {
				return rand.nextBoolean()
						? t
						: u;
			}
		};
	}
	
	/**
	 * CrossOverOperator that finds the single longest common subsequence of the two parent genomes (with the given amount of mismatches allowed).
	 * The resulting Sequence after cross-over contains this common sequence plus the sequence to the left taken from a random parent, and the sequence to the right taken from a random parent.
	 * Note: 1/2 probability that this operator simply returns one of the parents.
	 */
	public static CrossOverOperator singleSynapse() {
		return new CrossOverOperator() {
			Random rand = new Random();
			@Override
			public Sequence apply(Sequence t, Sequence u) {
				Pair<Sequence.Subsequence, Sequence.Subsequence> lcss = Sequence.LCSS(t, u);
				Sequence left = rand.nextBoolean()
						? t.get(0, lcss.left.start)
						: u.get(0, lcss.right.start);
				Sequence middle = t.get(lcss.left.start, lcss.left.start+lcss.left.length());
				Sequence right = rand.nextBoolean()
						? t.get(lcss.left.start + lcss.left.length(), t.length())
						: u.get(lcss.right.start + lcss.left.length(), u.length());
				return Sequence.cat(left, middle, right);
			}
		};
	}
	
	/**
	 * The Synapsing Variable Length Crossover algorithm as described in Hutt & Warwick 2007.
	 * @param 	tresholdLength
	 * 			The minimum length of a synapse.
	 */
	public static CrossOverOperator SVLC(int tresholdLength) {
		return new CrossOverOperator() {
			@Override
			public Sequence apply(Sequence t, Sequence u) {
				return Sequence.SVLC(t, u, tresholdLength).left;
			}
		};
	}

}
