package experiments;

import java.util.ArrayList;
import java.util.List;

import general.Distribution;
import organisms.Text;
import physics.CoordinateSystem;
import physics.TrivialCoordinateSystem;
import physics.World;
import settings.CrossOverOperator;
import settings.FitnessFunction;
import settings.LocationFilter;
import settings.MutationOperator;
import settings.SelectionRule;
import settings.Settings;

public class _5_TextEvolution {
	
	public static void main(String[] args) {
		
		// Initialization settings
		int worldSize = 400;
		int initialGenomeLength = 1500;
		int populationSize = 200;
		int experimentLength = 2000000;
		
		CoordinateSystem coordinateSystem = new TrivialCoordinateSystem(worldSize);
		Distribution<Integer> initialLengthDistribution = Distribution.constant(initialGenomeLength);
	
		// Evolution settings
		String targetText = "Shall I compare thee to a summer's day?";
		double wordGapPenalty = 5.;
		double letterGapPenalty = 1.;
		int minimumSynapseLength = 8*3;
		double substitutionProbability = .0002;
		double indelProbability = 0.0001;
		double splitProportion = 0.2;
		double sexyProportion = 0.01;
		
		MutationOperator mutationOperator = MutationOperator.combine(
				MutationOperator.poisson_sub(substitutionProbability),
				MutationOperator.poisson_indel(indelProbability,0,0,0,0,0,0,indelProbability,0,0,0,0,0,0,0,indelProbability)
				);
		CrossOverOperator crossoverOperator = CrossOverOperator.SVLC(minimumSynapseLength);
		Text.Factory factory = new Text.Factory(initialLengthDistribution, mutationOperator, crossoverOperator);
		FitnessFunction<Text> fitnessFunction =
				FitnessFunction.power(
				FitnessFunction.cutoff(
						Text.textComparison(targetText, wordGapPenalty, letterGapPenalty),
				0.,1.,true),
				10.);
		
		// Initialize
		Settings<Text> settings = new Settings<>();
		settings.setCoordinateSystem(coordinateSystem);
		settings.setOrganismFactory(factory);
		settings.setSpawnRule(LocationFilter.empty(), 1);
		settings.setSplittingRule(SelectionRule.selectProportion(fitnessFunction, splitProportion));
//		settings.setSexyRule(PairSelectionRule.any(SelectionRule.selectProportion(fitnessFunction,sexyProportion)));
		settings.setDeathRule(SelectionRule.cullToNumber(populationSize));
		
		List<Text> initialPopulation = new ArrayList<>();
		while (initialPopulation.size() < populationSize) {
			Text t = factory.random();
			if (fitnessFunction.apply(t) > 10.) {
				System.err.println(fitnessFunction.apply(t) + " <-- " + t);
				initialPopulation.add(t);
			}
		}
		
		World<Text> w = new World<>(settings);
		for (Text t : initialPopulation) w.spawnOrganism(t, w.getRandomEmptyLocation());

		// Evolve
		for (int i=0; i<experimentLength; i++) {
			w.tick();
//			System.err.println(Text.TextComparison(targetText, gapPenalty).apply(w).values());
			printInfo(w, fitnessFunction);
		}
		
	}
	
	public static void printInfo(World<Text> w, FitnessFunction<Text> f) {
		System.out.println(
				 w.getTime()
		+ ", " + w.getPopulation().size()
		+ ", " + f.totalFitness(w)
		+ ", " + w.getRandomOrganism()
		);
	}

}
