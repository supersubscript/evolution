package general;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;

/**
 * Abstract class that represents a probability distribution.
 * Provides several methods that return implementations.
 * @param N
 * 		The domain of this Distribution.
 */
@FunctionalInterface
public interface Distribution<N> extends Supplier<N> {
	
	/**
	 * Draw a random element from this this distribution.
	 */
	@Override
	public default N get() {
		return this.getRandom();
	}

	public N getRandom();
	
	/**
	 * Draw amt random elements from this distribution (with replacement).
	 */
	public default Collection<N> getRandom(int amt) {
		Collection<N> result = new ArrayList<>();
		for (int i=0; i<amt; i++) result.add(this.getRandom());
		return result;
	}
	
	/**
	 * Draw amt random elements from this distribution (without replacement)
	 * NOTE: default implementation is a simple while loop and may be inefficient or even result in an endless loop for Distributions with finite domains
	 * @throws IllegalArgumentException
	 * 			If the domain of this Distribution is finite and amt is greater than its cardinality.
	 */
	public default Collection<N> getRandomWithoutReplacement(int amt) {
		Collection<N> result = new HashSet<>();
		while (result.size() < amt) result.add(this.get());
		return result;
	}
	
	public static Distribution<Integer> binomial(int n, double p) {
		return new Distribution<Integer>() {
			Random rand = new Random();
			@Override
			public Integer getRandom() {
			  double log_q = Math.log(1.0 - p);
			   int x = 0;
			   double sum = 0;
			   for(;;) {
			      sum += Math.log(rand.nextDouble()) / (n - x);
			      if(sum < log_q) {
			         return x;
			      }
			      x++;
			   }
			}
		};
	}
	
	public static <N> Distribution<N> constant(N n) {
		return new Distribution<N>() {
			@Override
			public N getRandom() {
				return n;
			}
			@Override
			public Collection<N> getRandomWithoutReplacement(int amt) {
				switch (amt) {
				case 0:
					return new HashSet<>();
				case 1:
					Set<N> result = new HashSet<>();
					result.add(n);
					return result;
				default:
					throw new IllegalArgumentException("Constant distribution cannot return more than one distinct element");
				}
			}
		};
	}
	
	public static Distribution<Double> normal(double mu, double sigma) {
		return new Distribution<Double>() {
			Random rand = new Random();
			@Override
			public Double getRandom() {
				return rand.nextGaussian()*sigma + mu;
			}
		};
	}
	
	public static Distribution<Integer> poisson(double lambda) {
		return new Distribution<Integer>() {
			Random rand = new Random();
			@Override
			public Integer getRandom() {
			  double L = Math.exp(-lambda);
			  double p = 1.0;
			  int k = 0;

			  do {
			    k++;
			    p *= rand.nextDouble();
			  } while (p > L);

			  return (int) (k - 1);
			}
		};
	}
	
	/**
	 * A uniform Double Distribution across the specified domain.
	 */
	public static Distribution<Double> uniform(double start, double end) {
		return new Distribution<Double>() {
			Random rand = new Random();
			@Override
			public Double getRandom() {
				return rand.nextDouble() * (end-start) + start;
			}
		};
	}
	
	/**
	 * A uniform Integer Distribution across the specified domain.
	 * Note: Repeatedly taking random numbers without replacement where amt is close to the domain size is inefficient; use generic uniform on a predefined set of integers instead.
	 */
	public static Distribution<Integer> uniform(int start, int end) {
		return new Distribution<Integer>() {
			Random rand = new Random();
			@Override
			public Integer getRandom() {
				return rand.nextInt(end-start) + start;
			}
			@Override
			public Collection<Integer> getRandomWithoutReplacement(int amt) {
				if (amt > 3*(end-start)) {
					return Distribution.super.getRandomWithoutReplacement(amt);
				} else {
					Set<Integer> domain = new HashSet<>();
					for (int i=start; i<end; i++) domain.add(i);
					return Distribution.uniform(domain).getRandomWithoutReplacement(amt);
				}
			}
		};
	}
	
	/**
	 * A uniform distribution over the elements of an array.
	 */
	public static <N> Distribution<N> uniform(N[] ns) {
		return new Distribution<N>() {
			Random rand = new Random();
			@Override
			public N getRandom() {
				return ns[rand.nextInt(ns.length)];
			}
			@Override
			public Collection<N> getRandomWithoutReplacement(int amt) {
				if (amt > 3*(ns.length)) {
					return Distribution.super.getRandomWithoutReplacement(amt);
				} else {
					List<N> domain = Arrays.asList(ns);
					Set<N> result = new HashSet<>();
					while (result.size() < amt) {
						int index = rand.nextInt(domain.size());
						result.add(domain.get(index));
						domain.remove(index);
					};
					return result;
				}
			}
		};
	}
	
	/**
	 * A uniform distribution over the elements of a Collection.
	 */
	public static <N> Distribution<N> uniform(Collection<N> ns) {
		return new Distribution<N>() {
			Random rand = new Random();
			List<N> targets = new ArrayList<N>(ns);
			@Override
			public N getRandom() {
				if (targets.size() == 0) return null;
				return targets.get(rand.nextInt(ns.size()));
			}
			@Override
			public Collection<N> getRandomWithoutReplacement(int amt) {
				if (amt > 3*(ns.size())) {
					return Distribution.super.getRandomWithoutReplacement(amt);
				} else {
					List<N> domain = new ArrayList<>(ns);
					Set<N> result = new HashSet<>();
					while (result.size() < amt) {
						int index = rand.nextInt(domain.size());
						result.add(domain.get(index));
						domain.remove(index);
					};
					return result;
				}
			}
		};
	}
	
	/**
	 * A finite distribution where the probability of an element being drawn is given by the provided weights map.
	 * Negative weights are regarded as zero.
	 * If all weights are zero, the result is a uniform distribution over all the key elements of the weights map.
	 * If the supplied map is empty or null, returns the constant distribution for null.
	 */
	public static <N> Distribution<N> weighted(Map<N,Double> weights) {
		if (weights == null || weights.isEmpty()) return Distribution.constant(null);
		
		double sum = 0;
		for (double d : weights.values()) sum += Math.max(0, d);
		if (sum == 0.) return Distribution.uniform(weights.keySet());
		
		NavigableMap<Double, N> partialWeights = new TreeMap<>();
		double cumsum = 0.;
		for (Map.Entry<N,Double> entry : weights.entrySet()) {
			if (entry.getValue() <= Double.MIN_VALUE) continue;
			cumsum += entry.getValue() / sum;
			partialWeights.put(cumsum, entry.getKey());
		}
		
		return new Distribution<N>() {
			Random rand = new Random();
			@Override
			public N getRandom() {
				return partialWeights.ceilingEntry(rand.nextDouble()).getValue();
			}
			// TODO override getRandomWithoutReplacement
			@Override
			public String toString() {
				return "Weighted distribution: "+weights;
			}
		};
	}
	
	/**
	 * A weighted probability distribution where the probability of each element x in ns corresponds to the value in weights with the same index.
	 */
	public static <N> Distribution<N> weighted(List<N> ns, List<Double> weights) {
		if (ns == null || weights == null || ns.isEmpty() || weights.isEmpty()) return constant(null);
		Map<N,Double> map = new HashMap<>();
		for (int i=0; i<Math.min(ns.size(), weights.size()); i++) map.put(ns.get(i), weights.get(i));
		return Distribution.weighted(map);
	}
	
	/**
	 * A weighted distribution where the probability of each element x in ns being drawn is f(x) 
	 */
	public static <N> Distribution<N> weighted(Set<N> ns, ToDoubleFunction<N> f) {
		if (ns == null || f == null || ns.isEmpty()) return constant(null);
		Map<N,Double> map = new HashMap<>();
		for (N n : ns) map.put(n, f.applyAsDouble(n));
		return Distribution.weighted(map);
	}
	
	/**
	 * A weighted distribution where the probability of each element x in ns being drawn is f(x) 
	 */	
	public static <N> Distribution<N> weighted(Set<N> ns, Function<N,Double> f) {
		if (ns == null || f == null || ns.isEmpty()) return constant(null);
		Map<N,Double> map = new HashMap<>();
		for (N n : ns) map.put(n, f.apply(n));
		return Distribution.weighted(map);
	}

}
