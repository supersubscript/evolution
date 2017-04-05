#ifndef DIGORG_BINOM_H
#define DIGORG_BINOM_H

#include <algorithm>

// Finds the p that gives a cumulative binomial distribution (p, n) from 0 to k
// equal to pvalue. If upperlimit is true, the distribution is instead taken
// from n-k to n.
double binom_p_limit(unsigned k, unsigned n, double pvalue, bool upperlimit);

// A pair of the upper and lower p limits from the above function.
std::pair<double, double> binom_p_limits(unsigned k, unsigned n, double pvalue);

// The above rounded to integers with a factor other_n
std::pair<unsigned, unsigned> binom_p_limits_int(unsigned k, unsigned n,
	double pvalue, int other_n);

#endif
