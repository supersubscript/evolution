package organisms;


/**
 * A trivial implementation of the Organism class.
 * Contains no internal state nor a genome of any kind.
 */
public class Trivial extends Organism<Trivial> {
	public Trivial() {
	}
	
	public static class Factory implements Organism.Factory<Trivial> {
		@Override
		public Trivial random() {
			return new Trivial();
		}
		@Override
		public Trivial copy(Trivial o) {
			return new Trivial();
		}
		@Override
		public Trivial split(Trivial o) {
			return new Trivial();
		}
		@Override
		public Trivial sex(Trivial mommy, Trivial daddy) {
			return new Trivial();
		}
	}

}
