package settings;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;

import organisms.Organism;
import physics.CoordinateSystem;

/**
 * Class representing the settings of a particular experimental run.
 * Several settings pertain to certain functionality, and are implemented in abstract classes.
 * These abstract classes work similar to functional interfaces, but are implemented as abstract classes because they need to be non-static to access the values of other settings, such as grid size.
 * Factory methods are included but custom rules can be defined by extending the relevant class.
 */
public class Settings<O extends Organism<O>> {	
	/**
	 * Random number generator.
	 */
	public final Random rand = new Random();
	
	private final Collection<SelectionRule<O>> deathRules = new ArrayList<>();
	
	/**
	 * Create a death rule.
	 * The death rule is defined as a SelectionRule; each tick, every Organism selected by the rule will die.
	 */
	public void setDeathRule(SelectionRule<O> rule) {
		this.deathRules.add(rule);
	}
	
	public Collection<SelectionRule<O>> getDeathRules() {
		return new ArrayList<>(this.deathRules);
	}
	
	private final Collection<SelectionRule<O>> splittingRules = new ArrayList<>();
	
	/**
	 * Create a splitting rule.
	 * The splitting rule is defined as a SelectionRule; each tick, every Organism selected by the rule will try to clone itself (with mutations).
	 * Multiple rules can be set.
	 */
	public void setSplittingRule(SelectionRule<O> rule) {
		if (rule == null) throw new IllegalArgumentException("Cannot set null as splitting rule");
		this.splittingRules.add(rule);
	}
	
	public Collection<SelectionRule<O>> getSplittingRules() {
		return this.splittingRules;
	}
	
	private final Collection<PairSelectionRule<O>> sexyRules = new ArrayList<>();
	
	/**
	 * Set the rules for sex in this iteration.
	 * In each tick, all pairs of Organisms that are selected according to these rules will try to reproduce sexually.
	 */
	public void setSexyRule(PairSelectionRule<O> rule) {
		if (rule == null) throw new IllegalArgumentException("Cannot set null as sexy rule");
		this.sexyRules.add(rule);
	}
	
	public Collection<PairSelectionRule<O>> getSexyRules() {
		return this.sexyRules;
	}

	
	private int spawnDistance = 1;
	private LocationFilter spawnFilter;
	
	/**
	 * Set the rules for spawning new Organisms.
	 * When a new Organism is created, it is placed in a random new cell in the neighborhood (max. distance given) of its parent (for asexual reproduction) or the intersection of its parents' neighborhoods (for sexual reproduction),
	 * that fulfills the given LocationFilter requirement.
	 */
	public void setSpawnRule(LocationFilter spawnFilter, int distance) {
		if (distance <= 0) throw new IllegalArgumentException("Spawn distance must be positive");
		if (spawnFilter == null) throw new IllegalArgumentException("Spawn filter cannot be null");
		this.spawnDistance = 1;
		this.spawnFilter = spawnFilter;
	}
	
	public LocationFilter getSpawnFilter() {
		if (this.spawnFilter ==  null) throw new IllegalStateException("Spawn filter has not been set");
		return this.spawnFilter;
	}
	
	public int getSpawnDistance() {
		if (this.spawnDistance == 0) throw new IllegalStateException("Spawn distance has not been set");
		return this.spawnDistance;
	}
	
	private Organism.Factory<O> organismFactory;
	
	/**
	 * Set the Factory used to generate Organisms in this run.
	 */
	public void setOrganismFactory(Organism.Factory<O> factory) {
		this.organismFactory = factory;
	}
	
	public Organism.Factory<O> getOrganismFactory() {
		if (this.organismFactory == null) throw new IllegalStateException("Organism factory has not been set");
		return this.organismFactory;
	}
	
	private CoordinateSystem coordinateSystem;
	
	public CoordinateSystem getCoordinateSystem() {
		if (this.coordinateSystem == null) throw new IllegalStateException("No coordinate system was set");
		return this.coordinateSystem;
	}
	
	public void setCoordinateSystem(CoordinateSystem coordinateSystem) {
		this.coordinateSystem = coordinateSystem;
	}
	
	private DaytimeRule<O> daytimeRuleRule = DaytimeRule.trivial();
	
	public void setDaytimeRuleRule(DaytimeRule<O> rule) {
		if (rule == null) throw new IllegalArgumentException("Cannot set null as daytime rule");
		this.daytimeRuleRule = rule;
	}
	
	public DaytimeRule<O> getDaytimeRuleRule() {
		return this.daytimeRuleRule;
	}
	
	public static void main(String[] args) {
		List<Double> testList = new ArrayList<>();
		testList.add(0.);
		testList.add(0.5);
		testList.add(1.);
		testList.add(-5.);
		testList.sort(null);
		System.out.println(testList);
	}
}
