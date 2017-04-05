package organisms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import general.Distribution;
import general.Sequence;
import settings.CrossOverOperator;
import settings.FitnessFunction;
import settings.MutationOperator;

/**
 * An Organism class whose phenotype is a human-readable text.
 * The binary genome is read as a text phenotype by identifying substrings between certain start and stop sequences (in-frame) with ASCII words.
 */
public class Text extends Organism<Text> {
	
	private static final Sequence spaceASCII = Sequence.valueOf("{");
	private static final Sequence endASCII = Sequence.valueOf("}");
	
	private Text(Sequence genome) {
		this.genome = genome;
		this.words = genome.startstopMatchCodons(spaceASCII, endASCII).stream().map(Sequence::toASCII).collect(Collectors.toList());;
	}
	
	private final Sequence genome;
	
	private final List<String> words;
	
	public List<String> words() {
		return new ArrayList<>(this.words);
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append('(');
		sb.append(this.words.size());
		sb.append('/');
		sb.append(this.genome.length());
		sb.append(") ");
		
		for (String w : this.words()) {
			sb.append(w.replaceAll("[^a-zA-Z0-9_]", "#"));
			sb.append(" - ");
		}
		return sb.toString();
	}
	
	public static class Factory implements Organism.Factory<Text> {
		
		private final Distribution<Integer> initialLength;
		private final MutationOperator mut;
		private final CrossOverOperator crossover;
		
		public Factory(Distribution<Integer> initialLength, MutationOperator mut, CrossOverOperator crossover) {
			this.initialLength = initialLength;
			this.mut = mut;
			this.crossover = crossover;
		}

		@Override
		public Text random() {
			Sequence seq = Sequence.random(this.initialLength.get());
			return new Text(seq);
		}

		@Override
		public Text copy(Text o) {
			return new Text(o.genome);
		}

		@Override
		public Text split(Text o) {
			Sequence seq = this.mut.mutate(o.genome);
			return new Text(seq);
		}

		@Override
		public Text sex(Text mommy, Text daddy) {
			Sequence seq = this.crossover.apply(mommy.genome, daddy.genome);
			return new Text(seq);
		}
		
	}
	
	/**
	 * FitnessFunction that results in higher fitness if the Text's words are more similar to the given text.
	 * That is, the words are extracted from the given Text, and a Needleman-Wunsch global alignment is used to compare the Text.words() word sequence to the provided target text.
	 * The dissimilarity measure used is the edit distance normalized by the length of the longest word.
	 * Gap penalty is linear.
	 */
	public static FitnessFunction<Text> textComparison(String target, double wordGapPenalty, double letterGapPenalty) {
		return new FitnessFunction<Text>() {
			
			List<String> targetText = Arrays.asList(target.split("[ \n\r]+"));
			
			@Override
			public Double apply(Text o) {
				
				List<String> thisText = o.words();
				double[][] scoreMatrix = new double[thisText.size()+1][targetText.size()+1];
				
				for (int i = 0; i < thisText.size()+1; i++) {
					for (int j = 0; j < targetText.size()+1; j++) {
						if (i==0&&j==0) continue;
						double ins = i==0 ? -Double.MAX_VALUE : scoreMatrix[i-1][j] - wordGapPenalty;
						double del = j==0 ? -Double.MAX_VALUE : scoreMatrix[i][j-1] - wordGapPenalty;
						double match = i==0||j==0 ? -Double.MAX_VALUE : scoreMatrix[i-1][j-1] + wordSimilarity(targetText.get(j-1), thisText.get(i-1));
						scoreMatrix[i][j] = Math.max(Math.max(ins, del), match);
					}
				}
				
				return scoreMatrix[thisText.size()][targetText.size()] + wordGapPenalty * targetText.size();
			}
			
			double wordSimilarity(String ref, String comp) {
				char[] refchars = ref.toCharArray();
				char[] compchars = comp.toCharArray();
				double[][] scoreMatrix = new double[refchars.length+1][compchars.length+1];
				
				for (int i = 0; i < refchars.length+1; i++)
					for (int j=0; j<compchars.length+1; j++) {
						if (i==0&&j==0) continue;
						double ins = i==0 ? -Double.MAX_VALUE : scoreMatrix[i-1][j] - letterGapPenalty;
						double del = j==0 ? -Double.MAX_VALUE : scoreMatrix[i][j-1] - letterGapPenalty;
						double match = i==0||j==0 ? -Double.MAX_VALUE : scoreMatrix[i-1][j-1] + letterSimilarity(refchars[i-1], compchars[j-1]);
						scoreMatrix[i][j] = Math.max(Math.max(ins, del), match);
					}
				
				return scoreMatrix[refchars.length][compchars.length] + letterGapPenalty * refchars.length;
			}
			
			double letterSimilarity(char a, char b) {
				if (a==b) return 1.;
				if (Character.toLowerCase(a) == Character.toLowerCase(b)) return 0.5;
				return 0.;
			}
			
		};
	}
	
	public static void main(String[] args) {
		Text a = new Text(Sequence.valueOf("{Shall}{I}{compare}{thee}{to}{a}{summer's}{day?*} "));
		Text b = new Text(Sequence.valueOf("{4rO} "));
		System.out.println(a);
		System.out.println(b);
		System.out.println(textComparison("Shall I compare thee to a summer's day", 5., 1.).apply(b));
	}

}
