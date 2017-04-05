package general;

/**
 * A class representing pairs (2-tuples) of objects.
 */
public class Pair<A, B> {
	public final A left;
	public final B right;
	
	public Pair(A left, B right) {
		this.left = left;
		this.right = right;
	}
	
	public Pair<B, A> swap() {
		return new Pair<>(this.right, this.left);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Pair)) return false;
		@SuppressWarnings("unchecked")
		Pair<A,B> casted = (Pair<A,B>)obj;
		return this.left.equals(casted.left) && this.right.equals(casted.right);
	}
	
	@Override
	public int hashCode() {
		int hash = 17;
		hash = hash * 31 + this.left.hashCode();
		hash = hash * 31 + this.right.hashCode();
		return hash;
//		return this.left.hashCode()*this.right.hashCode() - this.right.hashCode();
	}
	
	@Override
	public String toString() {
		return "(" + this.left + ", " + this.right + ")";
	}
	
	
	
}
