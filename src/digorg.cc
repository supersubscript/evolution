#include <sstream>
#include <vector>
#include <string>
#include <iostream>

#include "common.h"
#include "population.h"
#include "recordlocker.h"

gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);

const long save_every = 1000;
const long log_every = 10;


void create_basedir(std::string basedir, size_t genome_length_init,
	std::string cost_type)
{
	organism::cost_type ct;
	if(cost_type == "ld")
		ct = organism::cost_clock_ld;
	else if(cost_type == "hf")
		ct = organism::cost_clock_hf;
	else
		throw "unknown cost type " + cost_type;

	organism::set_cost_type(ct);
	organism::set_genome_length_init(genome_length_init);

//	gene::set_promoter(pstart, pend);
//	gene::set_enable_repressors(repressors);

	if(mkdir(basedir.c_str(), 0777))
		throw "error on mkdir " + basedir;
	organism::save_globals(basedir);
}

void evolve_many(std::string basedir, size_t ensemble_size,
	size_t pop_size, size_t generations)
{
	organism::load_globals(basedir);

	recordlocker locker(basedir, ensemble_size);
	while(auto r = locker.lock_record(generations))
	{
		std::ostringstream dir;
		dir << basedir << "/" << r;
		population evolver(dir.str());
		if(!evolver.load())
		{
			evolver.resize(pop_size);
			std::cout << "Starting work on unit " << r << std::endl;
		}
		else
			std::cout << "Loaded unit " << r << std::endl;

		mkdir(dir.str().c_str(), 0777);
		std::ofstream log(dir.str() + "/evolution.log", std::ios::app);
		std::ofstream savelog(dir.str() + "/best.save", std::ios::app);
		std::ofstream cross_log(dir.str() + "/crossover.log", std::ios::app); 
		std::ofstream mut_log(dir.str() + "/mutation.log", std::ios::app); 
	   
      evolver.evolution(-generations, save_every, log_every, &log, &savelog, &cross_log, &mut_log);
		evolver.save();
		locker.unlock_record(evolver.generation());
		//evolver.print_best_timecourse();
	}
	std::cout << "All units done/started" << std::endl;
}


void anneal_many(std::string basedir, size_t ensemble_size,
					  size_t iterations, double t_start, double t_stop)
{
	organism::load_globals(basedir);

	recordlocker locker(basedir, ensemble_size);
	while(auto r = locker.lock_record())
	{
		std::ostringstream dir;
		dir << basedir << "/" << r;
		population evolver(dir.str());
		if(!evolver.load())
		{
			evolver.resize(2);  //anneal only needs two idividuals
			std::cout << "Starting work on unit " << r << std::endl;
		}
		else
			std::cout << "Loaded unit " << r << std::endl;

		mkdir(dir.str().c_str(), 0777);
		dir << "/anneal.log";
		std::ofstream log(dir.str(), std::ios::app);

		evolver.anneal(iterations, t_start, t_stop, &log, save_every, log_every);
		evolver.save();
		locker.unlock_record(evolver.generation());
		//evolver.print_best_timecourse();
	}
	std::cout << "All units done/started" << std::endl;
}

int printusage(char const* const* argv)
{
	std::string a = std::string(argv[0]);
	std::cerr <<
		"usage: " << a << " <basedir> init <initial genome size> ld/hf\n"
		"    or " << a << " <basedir> evolve_many <ensemblesize> <popsize> <generations>\n"
		"    or " << a << " <basedir> anneal_many <ensemblesize> <iterations> <temp_start> <t_stop>\n"
		"    or " << a << " <loaddir> load [optional: savedir]\n"
		"    or " << a << " <loaddir> eval\n"
		"    or " << a << " <loaddir> evolve_hist <prune_limit(-1/0/pos)>\n"
		"    or " << a << " <loaddir> randomnetworks <iterations> <temp_start> <t_stop> <savedir> [multiplier]\n";
	return 1;
}

int main(int argc, const char **argv)
{
	gsl_rng_set(rng, randrandseed());
//	gsl_rng_set(rng, 42);

	if(argc < 3)
		return printusage(argv);

	// directory to save/load from
	std::string modeldir(argv[1]);
	// Extract/remember the directory we're running from
	sourcedir(argv[0]);

	std::string mode(argv[2]);
	int args;      // correct length of argc for
	int optargs = 0;   // allow some optional arguments

	if(mode == "load")
	{
		args = 3;
		optargs = 1;
	}
	else if(mode == "randomnetworks")
	{
		args = 7;
		optargs = 1;
	}
	else if(mode == "evolve_many")
		args = 6;
	else if(mode == "anneal_many")
		args = 7;
	else if(mode == "init")
		args = 5;
	else if(mode == "evolve_hist")
		args = 4;
	else if(mode == "eval")
		args = 3;
	else
		return printusage(argv);

	if(argc < args || argc > args + optargs)
		return printusage(argv);

	try
	{
		if(mode == "load")
		{
			organism::load_globals(modeldir + "/..");
			population evolver(modeldir);
			if(!evolver.load())
				throw "failed to load from " + modeldir;

			// set new save path?
			if(argc > 3)
			{
				evolver.setdirectory(argv[3]);
				evolver.save();
			}
		}
		else if(mode == "randomnetworks")
		{
			int iterations = touint(argv[3]);
			double t_start = todouble(argv[4]);
			double t_stop  = todouble(argv[5]);

			std::string savedir = argv[6];
			mkdir(savedir.c_str(), 0777);
			std::cout << "savedir\t" << savedir << std::endl;

			int multiplier = 1; // how many nets to evolve per loaded net
			if(argc > 7)
				multiplier = touint(argv[7]);

			// Check if loadpath is a directory with many networks;
			// otherwise assume it's a single one.
			std::vector<std::string> dirnames;
			std::ifstream infile(std::string(modeldir + "/records"));
			if(infile)
			{
				// find which networks (populations) are done
				std::set<unsigned> recs;
				std::string a;
				while(infile >> a)
				{
					recs.insert(touint(a));
					infile >> a;
				}
				for(auto r : recs)
				{
					std::ostringstream tmp;
					tmp << modeldir << "/" << r;
					dirnames.push_back(tmp.str());
				}
				organism::load_globals(modeldir);
			}
			else
			{
				dirnames.push_back(modeldir);
				organism::load_globals(modeldir + "/..");
			}
			organism::save_globals(savedir);

			recordlocker locker(savedir, dirnames.size() * multiplier);

			while(auto r = locker.lock_record(1))
			{
				auto dir = dirnames[(r - 1) / multiplier];
				std::ostringstream subdir;
				subdir << savedir << "/" << r;

				std::cout << "doing:\t" << dir << " -> " << subdir.str() << std::endl;
				// the loaded population
				population orig_pop(dir);
				orig_pop.load();

				// find best organism
				double best_cost;
				organism &best = orig_pop.findbest(best_cost);
				std::cout << "best_organism:\t" << best_cost << std::endl;

				// the random network we want to evolve.
				int popsize = 2;
				population evolver(subdir.str(), popsize);

				mkdir(subdir.str().c_str(), 0777);
				subdir << "/anneal.log";
				std::ofstream log(subdir.str(), std::ios::app);

				evolver.anneal_null(best, iterations, t_start, t_stop, &log, save_every, log_every);
				evolver.save();
				locker.unlock_record(1);
			}
		}
		else if(mode == "evolve_hist")
		{
			double prune = todouble(argv[3]);
			organism::load_globals(modeldir);

			std::map<recordlocker::record_t, recordlocker::value_t> records;
			std::string infname = modeldir + "/records";
			std::ifstream infile(infname);
			if(!infile)
				throw "failed to open " + infname;
			// find which networks (populations) are done
			std::string a;
			while(getline(infile, a))
			{
				std::istringstream ss(a);
				recordlocker::record_t r;
				recordlocker::value_t v;
				ss >> r >> v;
				if(!ss || !ss.eof())
					throw "failed to read records from " + infname;
				records[r] = v;
			}
			if(records.empty())
				throw std::string("no records found");

			auto max_record = records.rbegin()->first;

			std::string evodir = modeldir + "/evohist_records";
			mkdir(evodir.c_str(), 0777);
			recordlocker locker(evodir, max_record);

			while(auto r = locker.lock_record())
			{
				// Skip unknown and unchanged populations
				if(records.find(r) == records.end() || records[r] == locker.get_value())
				{
					locker.unlock_record();
					continue;
				}

				std::ostringstream dir;
				dir << modeldir << "/" << r;
				std::cout << "doing:\t" << dir.str() << std::endl;

				population pop(dir.str());

				std::ostringstream fname;
				if(prune < 0)
					fname << dir.str() << "/best_unpruned.log";
				else
					fname << dir.str() << "/best_pruned_" << prune << ".log";

				std::ofstream log(fname.str());
				pop.loadbest(log, prune);

				locker.unlock_record(records[r]);
			}
			std::cout << "All units done/started" << std::endl;
		}
		else if(mode == "eval")
		{
			// todo: check that modeldir exists.

			organism::load_globals(modeldir + "/..");
			population evolver(modeldir);
			if(!evolver.load())
				throw "failed to load from " + modeldir;

			evolver.eval(std::cout);
		}
		else if(mode == "init")
		{
			create_basedir(modeldir, touint(argv[3]), argv[4]);
		}
		else if(mode == "evolve_many")
		{
			unsigned ensemble_size = touint(argv[3]);
			unsigned pop_size = touint(argv[4]);
			unsigned generations = touint(argv[5]);

			evolve_many(modeldir, ensemble_size, pop_size, generations);
		}
		else if(mode == "anneal_many")
		{
			unsigned ensemble_size = touint(argv[3]);
			unsigned iterations = touint(argv[4]);
			double t_start = todouble(argv[5]);
			double t_stop  = todouble(argv[6]);

			anneal_many(modeldir, ensemble_size, iterations, t_start, t_stop);
		}
		else
			throw std::string("mode error in main");
	}
	catch(std::string &s)
	{
		std::cerr << s << std::endl;
		return 1;
	}

	return 0;
}
