#include <gsl/gsl_cdf.h>
#include <iostream>
#include <cmath>
#include "binom.h"

double binom_p_limit(unsigned k, unsigned n, double pvalue, bool upperlimit)
{
	if(!n)
		return 0;
	double epsilon = 1e-9;
	double p;
	double plower = 0, pupper = 1;
	for(int i = 0; i < 100; ++i)
	{
		p = .5 * (pupper + plower);
		double pv = upperlimit ? gsl_cdf_binomial_P(k, p, n) : gsl_cdf_binomial_Q(k, p, n);
		if(std::abs(pv - pvalue) < epsilon * pvalue)
			break;
		if((pv > pvalue) == upperlimit)
			plower = p;
		else
			pupper = p;
	}
	return p;
}

std::pair<double, double> binom_p_limits(unsigned k, unsigned n, double pvalue)
{
	return std::make_pair(binom_p_limit(k, n, pvalue, false),
		binom_p_limit(k, n, pvalue, true));
}

std::pair<unsigned, unsigned> binom_p_limits_int(unsigned k, unsigned n,
	double pvalue, int other_n)
{
	return std::make_pair(
		(int)std::ceil(binom_p_limit(k, n, pvalue, false) * other_n),
		(int)std::floor(binom_p_limit(k, n, pvalue, true) * other_n));
}

/*
void print(unsigned k, unsigned n)
{
	auto lim = binom_p_limits(k, n, .01 * .5);
	std::cout << k << "\t" << n << "\t" << lim.first << "\t" <<
		lim.second << "\n";
}

int main()
{
	std::cout << "P = 0.01\n";
	std::cout << "#n\tk\tlower\tupper\n";

	print(107, 249);
	print(59, 249);
	print(0, 249);
	print(1, 249);
	print(67, 249);
	print(0, 249);
	print(1, 249);
	print(14, 249);

	print(281, 4536);
	print(523, 4536);
	print(0, 4536);
	print(1178, 4536);
	print(1127, 4536);
	print(0, 4536);
	print(574, 4536);
	print(853, 4536);
}
*/