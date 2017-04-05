package general;
import java.util.ArrayList;
import java.util.List;

/**
 * Abstract class that efficiently calculates and remembers primes.
 */
public abstract class Prime {
	
	private static final List<Integer> primes = new ArrayList<Integer>(); 
	
	public static int get(int i) {
		if (primes.isEmpty()) { primes.add(2); primes.add(3); };
		while (primes.size() <= i) {
			primes.add(nextPrime(primes.get(primes.size()-1)));
		}
		return primes.get(i);
	}
	
	public static int nextPrime(int p) {
		if (p<2) return 2;
		if (primes.contains(p) && primes.lastIndexOf(p) < primes.size()-1) return primes.get(primes.lastIndexOf(p)+1);
		int q = p+1;
		while (!primeQ(q)) q++;
		return q;
	}
	
	public static boolean primeQ(int p) {
		if (p <= primes.get(primes.size()-1)) return primes.contains(p);
		int i=0;
		while (get(i) <= Math.sqrt(p)) {
			if (p % get(i) == 0) return false;
			i++;
		}
		return true;
	}

}
