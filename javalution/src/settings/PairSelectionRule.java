package settings;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import general.Pair;
import organisms.Organism;
import physics.World;

/**
 * Strategy interface for selecting pairs of Organisms, such as when deciding which pairs are to reproduce sexually.
 */
@FunctionalInterface
public interface PairSelectionRule<O extends Organism<O>> extends Function<World<O>, Collection<Pair<O,O>>> {
	
	/**
	 * Alias of {@link #select(World)}.
	 */
	@Override
	default Collection<Pair<O, O>> apply(World<O> t) {
		return this.select(t);
	}

	/**
	 * Which Organisms to select, out of all the Organisms in the World.
	 */
	public Collection<Pair<O,O>> select(World<O> w);
	
	public static <O extends Organism<O>> PairSelectionRule<O> any(SelectionRule<O> rule) {
		return new PairSelectionRule<O>() {
			@Override
			public Collection<Pair<O, O>> select(World<O> w) {
				List<Pair<O,O>> result = new ArrayList<Pair<O,O>>();
				List<O> sel = new ArrayList<>(rule.select(w));
				Collections.shuffle(sel);
				for (int i=0; i<sel.size()-(sel.size()%2); i+=2) {
					result.add(new Pair<>(sel.get(i), sel.get(i+1)));
				}
				return result;
			}
		};
	}
	
	/**
	 * Selects any pair of Organisms that are selected by the supplied rule and also in each other's neighborhood.
	 * Once a pair has been selected, it is removed from the searching procedure, so that any Organism can occur in a selected pair according to this rule at most as many times at it is selected by the supplied SelectionRule.
	 */
	public static <O extends Organism<O>> PairSelectionRule<O> neighborhood(SelectionRule<O> rule, int distance) {
		return new PairSelectionRule<O>() {
			@Override
			public Collection<Pair<O, O>> select(World<O> w) {
				List<Pair<O,O>> result = new ArrayList<>();
				List<O> sel = new ArrayList<>(rule.select(w));
				Collections.shuffle(sel);
				iteratingSel:
				while (!sel.isEmpty()) {
					O first = sel.get(0);
					List<O> neighbors = new ArrayList<>(w.getNeighbors(first, distance));
					Collections.shuffle(neighbors);
					for (O second : neighbors) {
						if (sel.contains(second)) {
							result.add(new Pair<>(first,second));
							sel.remove(first);
							sel.remove(second);
							continue iteratingSel;
						}
					}
					sel.remove(first);
					continue iteratingSel;
				}
				return result;
			}
		};
	}
}