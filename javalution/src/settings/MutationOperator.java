package settings;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import general.Distribution;
import general.Sequence;

/**
 * This strategy interface that defines how a Sequence can be mutated.
*/
@FunctionalInterface
public interface MutationOperator {
	
	public Sequence mutate(Sequence seq);
	
	/**
	 * Trivial implementation. Does nothing.
	 */
	public static MutationOperator trivial() {
		return new MutationOperator() {
			@Override
			public Sequence mutate(Sequence seq) {
				return seq.clone();
			}
		};
	}
	
	/**
	 * Combine several mutations.
	 */
	public static MutationOperator combine(MutationOperator... parts) {
		return new MutationOperator() {
			@Override
			public Sequence mutate(Sequence seq) {
				Sequence result = seq;
				for (MutationOperator sub : parts) {
					result = sub.mutate(result);
				}
				return result;
			}
		};
	}
	
	/**
	 * A flip mutator that occurs randomly across the genome with a specified probability.
	 */
	public static MutationOperator poisson_sub(double prob) {
		return new MutationOperator() {
			@Override
			public Sequence mutate(Sequence seq) {
				return seq.xor(Sequence.random(seq.length(), prob));
			}
		};
	}
	
	/**
	 * A deletion mutator that occurs randomly across the genome.
	 * The probabilities are defined by a sequence; the first defines the probability of deleting 1 bit, the second defines the probability of deleting a run of two bits, etc.
	 * The probabilities are divided by their index to avoid bias towards deleting long runs (i.e. del(0, 0, 0.5) means that half of the bits will be deleted in runs of three; it does not mean that half of the indexes are the starting site of a deletion of length 3)
	 * When multiple random deletions overlap, the rightmost is shifted (i.e. if there is a deletion of length 2 at index 1,2 and another at index 2,3 then the result will be the deletion of the bits at index 1,2,3,4.)
	 * Deletion runs that would go over the length of the sequence are cut (i.e. if the sequence has length 10 and a deletion of length 3 is randomly chosen to start at index 8, then only indexes 8,9 will be deleted)  
	 */
	public static MutationOperator poisson_del(double... probs) {
		return new MutationOperator() {
			Random rand = new Random();
			@Override
			public Sequence mutate(Sequence seq) {
				List<Integer> remove = new ArrayList<Integer>();
				for (int i=0; i<probs.length; i++) {
					int amt = Distribution.binomial(seq.length()/(i+1), probs[i]).getRandom();
					
					for (int j=0; j<amt; j++) {
						int index = this.rand.nextInt(seq.length());
						int added = 0;
						while (added <= i && index < seq.length()) {
							if (!remove.contains(index)) {
								remove.add(index);
								added++;
							}
							index++;
						}
					}
				}
				Collections.sort(remove);

				boolean[] result = new boolean[seq.length() - remove.size()];
				// i is the index on the new, j the index on the old sequence
				int i=0; int j=0;
				while (i<result.length) {
					while (!remove.isEmpty() && remove.get(0) == j) {
						remove.remove(0);
						j++;
					}
					result[i] = seq.get(j);
					i++; j++;
				}
				return Sequence.create(result);
			}
		};
	}
	
	/**
	 * A duplication mutator that occurs randomly across the genome.
	 * The probabilities are defined by a sequence; the first defines the probability of duplicating 1 bit, the second defines the probability of duplicating a run of two bits, etc.
	 * The probabilities are divided by their index to avoid bias towards duplicating long runs. 
	 */
	public static MutationOperator poisson_in(double... probs) {
		return new MutationOperator() {
			Random rand = new Random();
			@Override
			public Sequence mutate(Sequence seq) {
				List<Integer> insert = new ArrayList<Integer>();
				for (int i=0; i<probs.length; i++) {
					int amt = Distribution.binomial(seq.length(), probs[i]/(i+1)).getRandom();
					
					for (int j=0; j<amt; j++) {
						int index = this.rand.nextInt(seq.length());
						for (int added=0; added<=i; added++) {
							insert.add(index);
						}
					}
				}
				Collections.sort(insert);

				boolean[] result = new boolean[seq.length() + insert.size()];
				// i is the index on the new, j the index on the old sequence
				int i=0; int j=0;
				while (i<result.length) {
					while (!insert.isEmpty() && insert.get(0) == j) {
						insert.remove(0);
						result[i] = rand.nextBoolean();
						i++;
					}
					result[i] = seq.get(j);
					i++; j++;
				}
				return Sequence.create(result);
			}
		};
	}

	/**
	 * Combination of insertions and deletion operators.
	 */
	public static MutationOperator poisson_indel(double... probs) {
		return MutationOperator.combine(
				MutationOperator.poisson_in(probs),
				MutationOperator.poisson_del(probs)
				);
	}

}
