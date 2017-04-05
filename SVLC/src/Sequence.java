import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public abstract class Sequence implements Comparable<Sequence> {
	
	private static final Random RAND = new Random();

	/**
	 * The Sequence length.
	 */
	public abstract int length();
	
	/**
	 * Get a particular bit from this Sequence.
	 * @return
	 */
	public abstract boolean get(int index);
	
	/**
	 * Clone this Sequence.
	 * Cloning a Sequence results in "flattening", i.e. a clone does not refer to any other Sequence objects.
	 */
	@Override
	public Sequence clone() {
		BitSet bits = new BitSet(this.length());
		for (int i=0; i<this.length(); i++) if (this.get(i)) bits.set(i);
		return Sequence.create(bits, this.length());
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("<");
		for (int i=0; i<this.length(); i++) sb.append(this.get(i) ? 1 : 0);
		sb.append(">");
		return sb.toString();
	}
	
	/**
	 * Compare two sequences.
	 * return 1 if this sequence is "larger" than argument (i.e. should be ordered first), -1 if smaller, 0 if equal.
	 * Comparison is typical dictionary ordering (with 1 > 0).
	 */
	@Override
	public int compareTo(Sequence o) {
		for (int i=0; i<this.length() && i<o.length(); i++) {
			if (this.get(i) && !o.get(i)) return 1;
			if (this.get(i) !=  o.get(i)) return -1;
		}
		if (this.length() < o.length()) return 1;
		if (this.length() > o.length()) return -1;
		return 0;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Sequence)) return false;
		return this.compareTo((Sequence)obj) == 0;
	}
	
	/**
	 * Create new Sequence from BitSet.
	 */
	public static Sequence create(BitSet bits, int length) {
		return new Standalonesequence(bits, length);
	}

	/**
	 * Create new Sequence from boolean array.
	 */
	public static Sequence create(boolean[] bits) {
		return new Standalonesequence(bits);
	}
	
	/**
	 * Create new Sequence from the concatenation of other Sequences.
	 */
	public static Sequence cat(Sequence... args) {
		int length = 0;
		BitSet bs = new BitSet();
		for (Sequence seq : args) {
			for (int i=0; i<seq.length(); i++) {
				if (seq.get(i)) bs.set(length+i);
			}
			length += seq.length();
		}
		return Sequence.create(bs, length);
	}

	/**
	 * Find the longest common subsequence of two sequences.
	 */
	public static Pair<Subsequence, Subsequence> LCSS(Sequence a, Sequence b) {
//		assert Sequence.LCSSdynamic(a, b).equals(Sequence.LCSSsuffixtree(a, b));
		return Sequence.LCSSsuffixtree(a,  b);
	}
	
	/**
	 * Find the longest common subsequence of two sequences.
	 * Implementation using dynamic programming.
	 * When comparing the first x characters in the UTF-8 binarized versions of Shakespear's Henry IV part 1 and 2, the least squares quadratic fit is:
	 * 		seconds = 0.017 - 0.008x + 0.011x^2
	 * (R2 = 0.9999)

	 */
	static Pair<Subsequence,Subsequence> LCSSdynamic(Sequence a, Sequence b) {
//		// Ensure a is longer than b
		if (a.length() < b.length()) return LCSSdynamic(b, a).swap();
		
		
		// initialize result
		int resultastart = 0;
		int resultbstart = 0;
		int resultlen = 0;
		
		// Slide sequence b across sequence a
		for (int astart = 0 ; astart < a.length() - resultlen; astart++) {
			int bstart = 0;
			int len = 0;
			while (astart+len < a.length() && bstart+len < b.length()) {
				len++;
				if (a.get(astart+len-1) == b.get(bstart+len-1)) {
					if (len > resultlen) {
						resultastart = astart;
						resultbstart = bstart;
						resultlen = len;
					}
				} else {
					bstart += len;
					len=0;
				}
			}
		}
		
		// Return result
		return new Pair<>(
				a.new Subsequence(resultastart, resultlen),
				b.new Subsequence(resultbstart, resultlen)
				);
	}
	
	/**
	 * Find the longest common subsequence of two sequences.
	 * Implementation using suffix tree.
	 * When comparing the first x characters in the UTF-8 binarized versions of Shakespear's Henry IV part 1 and 2, the least squares quadratic fit is:
	 * 		seconds = 0.008 + 0.004x + 0.00006x^2
	 * (R2 = 0.997)
	 */
	static Pair<Subsequence, Subsequence> LCSSsuffixtree(Sequence a, Sequence b) {

		// Generate and sort an array of "suffix" subsequences
		List<Subsequence> suffixarray = new ArrayList<>();
		for (int i=0; i<a.length(); i++) {
			suffixarray.add(a.new Subsequence(i, a.length()-i));
		}
		for (int i=0; i<b.length(); i++) {
			suffixarray.add(b.new Subsequence(i, b.length()-i));
		}
		Collections.sort(suffixarray);
		
		// Now the two substrings starting at the LCSS must be next to each other in the list
		// Find the two neighboring substrings with the highest startSimilarity
		int maxsim = 0;
		Subsequence asuffix = a.new Subsequence(0,0);
		Subsequence bsuffix = b.new Subsequence(0,0);
		for (int i=0; i<suffixarray.size()-1; i++) {
			if (suffixarray.get(i).supersequence() == suffixarray.get(i+1).supersequence()) continue;
			int sim = Sequence.startSimilarity(suffixarray.get(i), suffixarray.get(i+1));
			if (sim > maxsim) {
				maxsim = sim;
				asuffix = suffixarray.get(i).supersequence() == a ? suffixarray.get(i) : suffixarray.get(i+1);
				bsuffix = suffixarray.get(i).supersequence() == b ? suffixarray.get(i) : suffixarray.get(i+1);
			}
		}
		
		// Return the two nearest subsequences.
		return new Pair<>(
				a.new Subsequence(asuffix.start, maxsim),
				b.new Subsequence(bsuffix.start, maxsim)
				);
	}
	
	/**
	 * Return the number of identical bits at the beginning of two sequences.
	 */
	public static int startSimilarity(Sequence a, Sequence b) {
		int minlen = Math.min(a.length(), b.length());
		for (int i=0; i<minlen; i++) {
			if (a.get(i) != b.get(i)) return i;
		}
		return minlen;
	}
	
	/**
	 * Perform a cross-over of two Sequences using the Synapsing Variable-Length Crossover algorithm.
	 * A treshold length must be given that specifies the minimal synapse size.
	 */
	public static Pair<Sequence,Sequence> SVLC(Sequence a, Sequence b, int tresholdLength) {
		Pair<Subsequence, Subsequence> lcss = Sequence.LCSS(a, b);
		
		if (lcss.left.length < tresholdLength) {
			return new Pair<>(a, b);
			
		} else {
			
			Sequence aleft = a.new Subsequence(0, lcss.left.start);
			Sequence bleft = b.new Subsequence(0, lcss.right.start);
			Sequence aright = a.new Subsequence(lcss.left.start + lcss.left.length, a.length() - (lcss.left.start + lcss.left.length));
			Sequence bright = b.new Subsequence(lcss.right.start + lcss.right.length, b.length() - (lcss.right.start + lcss.right.length));
			
			Pair<Sequence, Sequence> leftSVLC = Sequence.SVLC(aleft, bleft, tresholdLength);
			Pair<Sequence, Sequence> rightSVLC = Sequence.SVLC(aright, bright, tresholdLength);
			boolean randleft = RAND.nextBoolean();
			boolean randright = RAND.nextBoolean();
			
			return new Pair<>(
					Sequence.cat(randleft ? leftSVLC.left : leftSVLC.right, lcss.left, randright ? rightSVLC.left : rightSVLC.right),
					Sequence.cat(randleft ? leftSVLC.right : leftSVLC.left, lcss.left, randright ? rightSVLC.right : rightSVLC.left)
					);
		}
	}
	
	/**
	 * 
	 */
	private static class Standalonesequence extends Sequence {
		/**
		 * Internal record of this Sequence's bit sequence.
		 */
		private final BitSet bits;
		
		/**
		 * Internal record of this Sequence's length (because BitSet only knows the index of the last 1)
		 */
		private final int length;
	
		/**
		 * Create new Sequence from BitSet.
		 */
		public Standalonesequence(BitSet bits, int length) {
			this.bits = bits;
			this.length = length;
		}

		/**
		 * Create new Sequence from boolean array.
		 */
		public Standalonesequence(boolean[] bits) {
			BitSet bitset = new BitSet(bits.length);
			for (int i = 0; i < bits.length ; i++) {
				if (bits[i]) bitset.set(i);
			}
			this.bits = bitset;
			this.length = bits.length;
		}
		
		@Override
		public boolean get(int index) {
			if (index >= this.length()) throw new IndexOutOfBoundsException();
			return this.bits.get(index);
		}
		
		@Override
		public int length() {
			return this.length;
		}
	}
	
	/**
	 * Convert a String to binary Sequence.
	 */
	public static Sequence valueOf(String s) {
		return Sequence.create(BitSet.valueOf(s.getBytes()), s.length()*8);
	}

	/**
	 * Class representing a subsequence of a sequence.
	 */
	public class Subsequence extends Sequence {
		
		public Subsequence(int start, int length) {
			if (start < 0) throw new IndexOutOfBoundsException("Subsequence must start at nonnegative index");
			if (length < 0) throw new IndexOutOfBoundsException("Subsequence must have nonnegative length");
			if (start+length > this.supersequence().length()) throw new IndexOutOfBoundsException("Subsequence cannot exceed supersequence");
			this.start = start;
			this.length = length;
		}
		
		public final int start;
		private final int length;
		
		public Sequence supersequence() {
			return Sequence.this;
		}
		
		@Override
		public boolean get(int index) {
			return Sequence.this.get(this.start + index);
		}
		
		@Override
		public int length() {
			return this.length;
		}
		
	}
	
	public static void main(String[] args) {
		Sequence a = Sequence.valueOf("The black cat lounged on the bed and meowed");
		Sequence b = Sequence.valueOf("The tortoiseshell cat sat on the chair and purred");
		System.out.println(a);
		System.out.println(b);
		
		System.out.println(Sequence.LCSSdynamic(a, b));
		System.out.println(Sequence.LCSSsuffixtree(a, b));
		System.out.println(Sequence.LCSS(a, b));
		
		
//		System.out.println(Sequence.SVLC(a, b, 8*3));
//		
//		System.out.println("The: " + Sequence.valueOf("The"));
//		System.out.println("cat: " + Sequence.valueOf("cat"));
//		System.out.println("on the: " + Sequence.valueOf("on the"));
//		System.out.println("and: " + Sequence.valueOf("and"));
	}

}
