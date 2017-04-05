package general;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A class representing a bijective association between specified keys and values.
 * This is essentially a simultaneous Map from A to B and from B to A.
 */
public class Bijection<A,B> {
	
	public Bijection() {
		this.mapA = new HashMap<A,B>();
		this.mapB = new HashMap<B,A>();
	}
	
	private final Map<A,B> mapA;
	private final Map<B,A> mapB;
	
	public void add(A a, B b) {
		if (this.hasA(a)) this.removeA(a);
		if (this.hasB(b)) this.removeB(b);
		this.mapA.put(a, b);
		this.mapB.put(b, a);
	}
	
	public void remove(A a, B b) {
		this.mapA.remove(a);
		this.mapB.remove(b);
	}
	
	public void removeA(A a) {
		if (this.hasA(a))
		this.remove(a, this.a2b(a));
	}
	
	public void removeB(B b) {
		if (this.hasB(b)) {
			this.remove(this.b2a(b), b);
		}
	}
	
	public boolean hasA(A a) {
		return this.mapA.containsKey(a);
	}
	
	public boolean hasB(B b) {
		return this.mapB.containsKey(b);
	}
	
	public B a2b(A a) {
		return this.mapA.get(a);
	}
	
	public A b2a(B b) {
		return this.mapB.get(b);
	}
	
	public List<B> a2b(List<A> as) {
		List<B> result = new ArrayList<>();
		for (A a : as) if (this.hasA(a)) result.add(this.a2b(a));
		return result;
	}
	
	public List<A> b2a(List<B> bs) {
		List<A> result = new ArrayList<>();
		for (B b : bs) if (this.hasB(b)) result.add(this.b2a(b));
		return result;
	}
	
	public Set<B> a2b(Set<A> as) {
		Set<B> result = new HashSet<>();
		for (A a : as) if (this.hasA(a)) result.add(this.a2b(a));
		return result;
	}
	
	public Set<A> b2a(Set<B> bs) {
		Set<A> result = new HashSet<>();
		for (B b : bs) if (this.hasB(b)) result.add(this.b2a(b));
		return result;
	}
	
	public Set<A> avalues() {
		return mapA.keySet();
	}
	
	public Set<B> bvalues() {
		return mapB.keySet();
	}
	
	@Override
	public String toString() {
		return this.mapA.toString();
	}

}
