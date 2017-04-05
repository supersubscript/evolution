package physics;

import java.util.HashSet;
import java.util.Set;

import general.Distribution;
import general.Pair;
import organisms.Organism;
import settings.PairSelectionRule;
import settings.SelectionRule;
import settings.Settings;

/**
 * A class representing a changing and evolving World.
 * A world is a Habitat that can be used as the basis for evolution. It has two added functionalities:
 * 		- reproduction: a World knows the rules for death and reproduction
 * 		- "time": the ability to change the internal state of all Organisms and kill and multiply
 * Any changes to this World should be made using the Action system to adhere to the rules of reproduction and time.
 * 
 * [In lieu of using an extensive interface and inheritance network, all extra functionality required for specific Organisms and experiments will simply be added to this class with a default implementation.]
 */
public class World<O extends Organism<O>>
	extends Habitat<O> {
	
	public World(Settings<O> settings) {
		super(settings.getCoordinateSystem());
		this.settings = settings;
	}
	
	public final Settings<O> settings;

	////// Reproduction
	
	/**
	 * Get all the possible locations for spawning a new Organism near the given coordinate.
	 */
	public Set<CoordinateSystem.Coordinate> getSpawnLocations(CoordinateSystem.Coordinate coordinate) {
		return CoordinateSystem.filter(
				coordinate.neighbors(this.settings.getSpawnDistance()),
				this.settings.getSpawnFilter().concretize(this)
				);
	}
	
	/**
	 * Get a random spawning location around the given coordinates.
	 * Return 0 if there is none.
	 */
	public CoordinateSystem.Coordinate getRandomSplitLocation(O o) {
		if (!this.contains(o)) return null;
		CoordinateSystem.Coordinate result = this.getLocation(o).randomNeighbor(
				this.settings.getSpawnDistance(),
				this.settings.getSpawnFilter().concretize(this));
		return result;
	}
	
	public CoordinateSystem.Coordinate getRandomHumpLocation(O o1, O o2) {
		if (!this.contains(o1)) return null;
		if (!this.contains(o2)) return null;
		Set<CoordinateSystem.Coordinate> neigh1 = this.getLocation(o1).neighbors(this.settings.getSpawnDistance());
		Set<CoordinateSystem.Coordinate> neigh2 = this.getLocation(o2).neighbors(this.settings.getSpawnDistance());
		Set<CoordinateSystem.Coordinate> intersection = neigh1; intersection.retainAll(neigh2);
		intersection = CoordinateSystem.filter(intersection, this.settings.getSpawnFilter().concretize(this));
		return Distribution.uniform(intersection).getRandom();
	}
	
	/**
	 * Spawn a random Organism in an empty location, if it exists.
	 */
	public void spawnRandom() {
		this.spawnOrganism(this.settings.getOrganismFactory().random(), this.getRandomEmptyLocation());
	}

	
	/////// Time
	
	private int time = 0;
	
	public int getTime() {
		return this.time;
	}
	
	/**
	 * Pass one turn.
	 * This invokes the tick() method for all Organisms, and invokes general death and replication rules.
	 */
	public void tick() {
		this.time++;
		
		// Make time pass for Organisms
		for (O o : this.getPopulation()) o.tick(this);
		
		// Kill dying Organisms
		for (SelectionRule<O> rule : this.settings.getDeathRules())
			for (O o : rule.select(this)) this.new Kill(o);
		
		// Split splitting Organisms
		for (SelectionRule<O> rule : this.settings.getSplittingRules())
			for (O o : rule.select(this)) this.new Split(o);
		
		// Reproduce sexy Organisms
		for (PairSelectionRule<O> rule : this.settings.getSexyRules())
			for (Pair<O,O> p : rule.select(this)) this.new Hump(p.left,p.right);
		
		// Execute Actions
		for (Action a : this.actions) a.execute();
		this.actions.clear();
	}
	
	public final Set<Action> actions = new HashSet<>(); 
	
	
	/**
	 * Abstract class representing an Action to be undertaken in the World.
	 * Actions provide an easy interface to make changes to the World while enforcing the rules of time and reproduction.
	 * In particular, Actions are undertaken at the end of each tick, ensuring, for example, that an Organism can reproduce and die in the same tick.
	 */
	public abstract class Action {
		protected Action() {
			World.this.actions.add(this);
		}
		protected abstract void execute();
	}
	
	/**
	 * Create a new Action that kills the target Organism at the end of the tick.
	 */
	public class Kill extends Action {
		public Kill(CoordinateSystem.Coordinate coordinates) {
			this(World.this.getOccupant(coordinates));
		}
		public Kill(O o) {
			this.organism = o;
		}
		private final O organism;
		
		@Override
		protected void execute() {
			World.this.kill(this.organism);
		}
	}
	
	/**
	 * Create a new Action that spawns a new Organism in the specified location at the end of the tick.
	 * Nothing happens if the target coordinate is already filled.
	 */
	public class Spawn extends Action {
		public Spawn(CoordinateSystem.Coordinate coordinates, O o) {
			this.coordinates = coordinates;
			this.organism = o;
		}
		private final CoordinateSystem.Coordinate coordinates;
		private final O organism;
		
		@Override
		protected void execute() {
			World.this.spawnOrganism(organism, coordinates);
		}
	}
	
	/**
	 * Create a new Action that splits the target Organism at the end of the tick.
	 * The spawn location is determined at the moment of creation for this Action.
	 * If it is filled at the moment of execution, no offspring will be created, but the Organism's Factory's split() method will still be called, which may change the inner state of the organism.
	 */
	public class Split extends Spawn {
		public Split(O o) {
			super(World.this.getRandomSplitLocation(o), World.this.settings.getOrganismFactory().split(o));
		}
	}
	
	/**
	 * Create a new Action that creates an offspring for the target Organisms.
	 * The spawn location is determined at the moment of creation for this Action.
	 * If it is filled at the moment of execution, no offspring will be created, but the Organism's Factory's sex() method will still be called, which may change the inner state of the organism.
	 */
	public class Hump extends Spawn {
		public Hump(O o1, O o2) {
			super(World.this.getRandomHumpLocation(o1, o2), World.this.settings.getOrganismFactory().sex(o1, o2));
		}
	}
	
	
	///// Light
	
	public double getLight(O o) {
		return this.settings.getDaytimeRuleRule().light(this, o);
	}

}
