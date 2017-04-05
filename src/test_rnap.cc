#include <vector>
#include <iostream>
#include <iterator>
#include "genes.h"

//#define DEBUG_GENES
#include "genes.cc"


gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);


int main(int argc, char **argv)
{
	if(argc != 4)
	{
		std::cerr <<
			"Tests if the optimized getrnapactivity gives the right result\n"
			"usage: " << argv[0] << " tfs iters randseed\n";
		return 1;
	}

	gsl_rng_set(rng, atoi(argv[3]));

	std::vector<gene> tfs(atoi(argv[1]));
	std::vector<double> levels(tfs.size(), .1);
//	gene::veclevelgetter getter(levels);
	gene reg;
	int iters = atoi(argv[2]);

	for(int i = 0; i < iters; ++i)
	{
		reg.randomize();
		for(auto &i: tfs)
			i.randomize();
		reg.updatesites(tfs);

		// double a = reg.getrnapactivity_safe(getter);
		// double b = reg.getrnapactivity(getter);
		// std::cout << 2 * std::abs(a - b) / (a + b) << "\n";
	}
}
