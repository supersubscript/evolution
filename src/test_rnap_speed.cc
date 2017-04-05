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
			"Benchmarks the (un)optimized getrnapactivity.\n"
			"usage: " << argv[0] << " tfs iters randseed\n";
		return 1;
	}

	gsl_rng_set(rng, atoi(argv[3]));

	std::vector<gene> tfs(atoi(argv[1]));
	for(auto &i: tfs)
		i.randomize();
	std::vector<double> levels(tfs.size(), .1);

	gene reg;
	reg.randomize();
	reg.updatesites(tfs);

	std::set<size_t> tfindexes;
	reg.getbindingtfs(tfindexes);
	std::cout << "Binding " << tfindexes.size() << " TFs\n";

	int iters = atoi(argv[2]);

	double a;
	getcputime();
//	for(int i = 0; i < iters; ++i)
//		a = reg.getrnapactivity(gene::veclevelgetter(levels));
	double t = getcputime();
	std::cout << "  optimized: " << a << "\t" << t << " sec" << std::endl;

//	for(int i = 0; i < iters; ++i)
//		a = reg.getrnapactivity_safe(gene::veclevelgetter(levels));
	t = getcputime();
	std::cout << "unoptimized: " << a << "\t" << t << " sec\n";
}
