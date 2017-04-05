package general;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Class representing a binary bit sequence.
 * This class is abstract to encapsulate the internal data storage, but concrete Sequences can be initialized using Sequence.create().
 */
public abstract class Sequence implements Comparable<Sequence> {
	
	private static final Random RAND = new Random();

	/**
	 * The Sequence length.
	 */
	public abstract int length();
	
	/**
	 * Get a particular bit from this Sequence.
	 */
	public abstract boolean get(int index);
	
	/**
	 * Get a particular Subsequence of this Sequence.
	 */
	public Subsequence get(int start, int stop) {
		return this.new Subsequence(start, stop-start);
	}
	
	/**
	 * Chop this Sequence into a List of Subsequences separated at the specified indices.
	 */
	public List<Subsequence> chop(List<Integer> indices) {
		List<Integer> ind = new ArrayList<>(indices);
		ind.add(this.length());
		Collections.sort(ind);
		List<Subsequence> result = new ArrayList<>();
		for (int i = 0; i < ind.size()-1; i++) {
			if (ind.get(i) <= 0 || ind.get(i) > this.length()) continue;
			result.add(this.get(ind.get(i), ind.get(i+1)));
		}
		return result;
	}
	
	/**
	 * @return this.chop(length, false)
	 */
	public List<Subsequence> chop(int length) {
		return this.chop(length, false);
	}
	
	/**
	 * Chop this Sequence into a list of Subsequences separated at every multiple of the given index.
	 * Overhang can be included or excluded.
	 */
	public List<Subsequence> chop(int length, boolean includeOverhang) {
		List<Subsequence> result = new ArrayList<>();
		int i = 0;
		while (i <= this.length()) {
			result.add(this.get(i, i+length));
			i += length;
		}
		if (includeOverhang && i < this.length()) result.add(this.get(i, this.length()));
		return result;
	}
	
	/**
	 * Split this Sequence by the occurrences of the specified sequence.
	 */
	public List<Subsequence> split(Sequence sep) {
		List<Subsequence> result = new ArrayList<>();
		int i=0;
		do {
			int seploc = this.find(sep, i);
			if (seploc == -1) break;
			result.add(this.get(0, seploc));
			i = seploc + sep.length();
		} while (i > 0);
		return result;
	}
	
	/**
	 * Find the location of a particular subsequence, only looking inside the specified subsequence index.
	 * The returned location is the index of the first bit of the first occurrence of the subsequence, of -1 if the subsequence is not present.
	 * Index starts from the start of this Sequence, not the start of the searching interval.
	 */
	public int find(Sequence needle, int start, int end) {
		int astartmax = Math.min(end, this.length()) - needle.length();
		astartloop:
		for (int astart = Math.max(0, start); astart < astartmax; astart++) {
			for (int j=0; j<needle.length(); j++) {
				if (this.get(astart+j) != needle.get(j)) continue astartloop;
			}
			return astart;
		}
		return -1;
	}
	
	/**
	 * @return this.find(needle, start, this.length())
	 */
	public int find(Sequence needle, int start) {
		return this.find(needle, start, this.length());
	}
	
	/**
	 * Find the location of a particular subsequence.
	 * The returned location is the index of the first bit of the first occurrence of the subsequence, of -1 if the subsequence is not present.
	 */
	public int find(Sequence needle) {
		return this.find(needle, 0, this.length());
	}
	
	/**
	 * Find the location of a particular subsequence, only searching in the specified reading frame.
	 * The returned location is the index of the first bit, or -1 if no match was found.
	 */
	public int findCodon(Sequence needle, int start, int frameSize) {
		if (start < 0) throw new IndexOutOfBoundsException("Start of codon search cannot be negative");
		int astartmax = this.length() - needle.length();
		astartloop:
		for (int astart = start; astart < astartmax; astart+=frameSize) {
			for (int j=0; j<needle.length(); j++) {
				if (this.get(astart+j) != needle.get(j)) continue astartloop;
			}
			return astart;
		}
		return -1;
	}
	
	/**
	 * Find the locations of all occurrences of a particular subsequence.
	 * This can be done in an overlapping or non-overlapping manner.
	 */
	public List<Integer> findAll(Sequence needle, boolean overlap) {
		List<Integer> result = new ArrayList<>();
		int startLooking = 0;
		int firstOcc = this.find(needle);
		while (firstOcc != -1) {
			result.add(firstOcc);
			startLooking = firstOcc + (overlap ? 1 : needle.length());
			firstOcc = this.find(needle, startLooking);
		}
		return result;
	}
	
	/**
	 * @return this.find(needle, true)
	 */
	public List<Integer> findAll(Sequence needle) {
		return this.findAll(needle, true);
	}
	
	/**
	 * Check if this Sequence contains a particular subsequence.
	 */
	public boolean contains(Sequence needle) {
		return this.find(needle) != -1;
	}
	
	/**
	 * Find matches for the pattern start...end.
	 * Returns all the bits between start and end sequences (exclusive).
	 * No overlap is allowed.
	 * Use this method for finding genes of freely variable length.
	 */
	public List<Subsequence> startstopMatch(Sequence start, Sequence stop) {
		List<Subsequence> result = new ArrayList<>();
		int startPos = -1;
		int endPos = -stop.length();
		while (true) {
			startPos = this.find(start, endPos+stop.length());
			if (startPos < 0) break;
			endPos = this.find(stop, startPos+start.length());
			if (endPos < 0) break;
			result.add(this.get(startPos+start.length(),endPos));
		}
		return result;
	}
	
	/**
	 * Find matches for the pattern start..end, where .. is a multiple of a fixed number of bits.
	 * Returns all the bits between start and end sequences (exclusive).
	 * No overlap is allowed.
	 * Use this method for finding genes composed of codons. 
	 */
	public List<Subsequence> startstopMatchCodons(Sequence start, Sequence stop) {
		if (start.length() != stop.length()) throw new IllegalArgumentException("Start and stop codons must be of equal length");
		List<Subsequence> result = new ArrayList<>();
		int startPos = -1;
		int endPos = -stop.length();
		while (true) {
			startPos = this.find(start, endPos+stop.length());
			if (startPos < 0) break;
			endPos = this.findCodon(stop, startPos+start.length(), start.length());
			if (endPos < 0) break;
			result.add(this.get(startPos+start.length(), endPos));
		}
		return result;
	}
	
	/**
	 * @return this.startstopMatch(start, stop)
	 */
	@Deprecated
	public List<Subsequence> genes(Sequence start, Sequence stop) {
		return this.startstopMatch(start, stop);
	}
	
	/**
	 * Get the complement of this sequence.
	 */
	public Sequence complement() {
		BitSet bits = new BitSet();
		for (int i=0; i<this.length(); i++) {
			bits.set(i, !this.get(i));
		}
		return Sequence.create(bits, this.length());
	}
	
	/**
	 * Get the reverse of this sequence.
	 */
	public Sequence reverse() {
		BitSet bits = new BitSet();
		for (int i=0; i<this.length(); i++) {
			bits.set(this.length()-i-1, i++);
		}
		return Sequence.create(bits, this.length());
	}
	
	/**
	 * Get the reverse complement of this sequence.
	 * @return this.reverse().complement()
	 */
	public Sequence reverseComplement() {
		BitSet bits = new BitSet();
		for (int i=0; i<this.length(); i++) {
			bits.set(this.length()-i-1, !this.get(i));
		}
		return Sequence.create(bits, this.length());
	}
	
	/**
	 * Return the sequence obtained by repeated boolean operation AND.
	 */
	public Sequence and(Sequence s) {
		if (this.length() != s.length()) throw new IllegalArgumentException("Can only do boolean comparison between Sequences of the same length");
		BitSet bits = new BitSet();
		for (int i=0; i<this.length(); i++) {
			bits.set(i, this.get(i) && s.get(i));
		}
		return Sequence.create(bits, this.length());
	}
	
	/**
	 * Return the sequence obtained by repeated boolean operation OR.
	 */
	public Sequence or(Sequence s) {
		if (this.length() != s.length()) throw new IllegalArgumentException("Can only do boolean comparison between Sequences of the same length");
		BitSet bits = new BitSet();
		for (int i=0; i<this.length(); i++) {
			bits.set(i, this.get(i) || s.get(i));
		}
		return Sequence.create(bits, this.length());
	}
	
	/**
	 * Return the sequence obtained by repeated boolean operation XOR.
	 */
	public Sequence xor(Sequence s) {
		if (this.length() != s.length()) throw new IllegalArgumentException("Can only do boolean comparison between Sequences of the same length");
		BitSet bits = new BitSet();
		for (int i=0; i<this.length(); i++) {
			bits.set(i, this.get(i) ^ s.get(i));
		}
		return Sequence.create(bits, this.length());
	}
	
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
	 * Create new random Sequence of given length.
	 */
	public static Sequence random(int length) {
		BitSet bits = new BitSet();
		for (int i=0; i<length; i++) {
			bits.set(i, RAND.nextBoolean());
		}
		return Sequence.create(bits, length);
	}
	
	/**
	 * Create a new random Sequence of given length, where 1 occurs with probability p.
	 */
	public static Sequence random(int length, double p) {
		BitSet bits = new BitSet();
		for (int i=0; i<length; i++) {
			bits.set(i, RAND.nextDouble() < p);
		}
		return Sequence.create(bits, length);
	}
	
	/**
	 * Create a new random Sequence of given length, where 1 occurs with probability p.
	 * Optimized for low values of p.
	 */
	public static Sequence randomLowP(int length, double p) {
		// TODO: Optimize using binomial distribution
		BitSet bits = new BitSet();
		for (int i=0; i<length; i++) {
			bits.set(i, RAND.nextDouble() < p);
		}
		return Sequence.create(bits, length);
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
		
		@Override
		public int hashCode() {
			return this.bits.hashCode() + this.length;
		}
	}
	
	/**
	 * Convert a String to binary Sequence using ASCII encoding.
	 */
	public static Sequence valueOf(String s) {
		return Sequence.create(BitSet.valueOf(s.getBytes()), s.length()*8);
	}
	
	/**
	 * Convert a Sequence to String using ASCII encoding.
	 */
	public String toASCII() {
		byte[] bytes = this.asBitSet().toByteArray();
		for (int i=0; i<bytes.length; i++) bytes[i] = (byte) (bytes[i] & ~0x80);
		return new String(bytes, StandardCharsets.US_ASCII);
	}
	
	/**
	 * Read this Sequence as a BitSet.
	 */
	public BitSet asBitSet() {
		BitSet result = new BitSet(this.length());
		for (int i=0; i<this.length(); i++) result.set(i, this.get(i));
		return result;
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
		
		@Override
		public int hashCode() {
			return this.start*this.supersequence().hashCode() + this.length;
		}
		
	}
	
	public static void main(String[] args) {
		for (byte b=-126; b<126; b++) {
			byte x = (byte) (b & ~0x80);
			System.out.println(b
					+ " -- " + x
					+ " -- " + String.format("%8s", Integer.toBinaryString(x & 0xFF)).replace(' ', '0')
					+ " -- " + new String(new byte[]{x}, StandardCharsets.US_ASCII)
					);
		}
	}

}
