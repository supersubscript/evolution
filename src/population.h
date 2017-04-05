#ifndef DIGORG_POPULATION_H

#include <algorithm>
#include <iterator>
#include <sys/stat.h>
#include <time.h>
#include "organism.h"


class population
{
public:
	// Creates a new population
	population(std::string dir) :
		_directory(dir), _generation(0), _cputime(0)
	{
	}

	population(std::string dir, size_t pop_size) :
		_directory(dir), _generation(0), _cputime(0)
	{
		resize(pop_size);
	}

	void resize(size_t pop_size)
	{
		if(pop_size < _organisms.size())
			_organisms.erase(_organisms.begin() + pop_size, _organisms.end());
		while(pop_size > _organisms.size())
			_organisms.push_back(organism(organism::get_genome_length_init()));
	}

	inline size_t generation() const
	{
		return _generation;
	}

	inline size_t cputime() const
	{
		return _cputime;
	}

	static void log_header(std::ostream &log, bool full)
	{
		log << "#generation\tcputime\tdiversity";
		organism::log_header(log, full);
	}

   static void log_header_anneal(std::ostream &log)
   {
      log << "#generation\tcost\tgenes\tlinks\ttemp" << std::endl;
	}

	void log_line(std::ostream &log, organism &best, bool full)
	{
		log << _generation << "\t" 
		<< _cputime << "\t" << get_pop_diversity(); 
		best.log_line(log, full);
	}

	void log_line(std::ostream &log, double temp, organism &best, bool full)
	{
		log << _generation << "\t" << temp;
		best.log_line(log, full);
	}

	// An evolutionary algorithm
	// if iterations is <0, use -iterations as the target generation number.
	void evolution(long iterations, long save_every, long log_every,
		std::ostream *log, std::ostream *savelog, std::ostream *cross_log, std::ostream *mut_log)
	{
		// Start counting CPU time
		getcputime();
		double time_between_saves = 15*60; // 15 minutes
		double next_save = _cputime + time_between_saves;
		
		simulator::log = nullptr;
		assert(!_organisms.empty());
		if(_organisms.size() < 3)
		{
			std::cerr << "evolution() requires popsize of at least 3.\n"
				"Adding random nets\n";
			resize(4);
		}

		if(log)
			log_header(*log, false);
      if(mut_log)
         *mut_log << "#generation\tbefore\tafter\treeval_after\tdeltagenes\n";
		std::vector<double> costs;
		for(size_t i = 0; i < _organisms.size(); ++i)
			costs.push_back(_organisms[i].getcost());

		size_t endgeneration;
		if(iterations < 0)
			endgeneration = -iterations;
		else
			endgeneration = _generation + iterations;

		for(; _generation < endgeneration; ++_generation)
		{
			// change one individual
			if(gsl_rng_uniform(rng) >= organism::get_gene_crossover_prob())
			{
				//replace worst individual (highest cost):
				int index_replace = std::distance(costs.begin(),
					std::max_element(costs.begin(), costs.end()));
			
		  		int index = choose_tournament(costs);
            _organisms[index_replace] = _organisms[index];
            _organisms[index_replace].mutate();
            double before = costs[index];
            double after = _organisms[index_replace].getcost();
            costs[index_replace] = after;
            if(std::abs(before-after) > std::numeric_limits<double>::epsilon() && mut_log)
            {
               *mut_log << _generation << "\t" << before;
               *mut_log << "\t" << after;

               costs[index_replace] = _organisms[index_replace].getcost(true);
               *mut_log << "\t" << costs[index_replace] << "\t" << 
                  int(_organisms[index_replace].size()) - int(_organisms[index].size()) 
                     <<  std::endl;
            }
         }
			else
			{
				std::vector<std::pair<double, size_t>> costpairs;
				
				for(size_t i = 0; i < costs.size(); ++i)
					costpairs.push_back(std::make_pair(costs[i], i));

				// sort to find the two worst individuals				
				std::sort(costpairs.begin(), costpairs.end());
				size_t index_replace1 = costpairs[costpairs.size()-1].second;
				size_t index_replace2 = costpairs[costpairs.size()-2].second;

				int index1 = choose_tournament(costs);
				int index2 = choose_tournament(costs);
		
				if (cross_log)
					*cross_log << _generation << "\t" << costs[index1] << "\t" << costs[index2];

				organism::crossover(_organisms[index1], _organisms[index2],
					_organisms[index_replace1], _organisms[index_replace2]);

				// update fitness
				costs[index_replace1] = _organisms[index_replace1].getcost();
				costs[index_replace2] = _organisms[index_replace2].getcost();

            // print data for crossover
				if(cross_log)
				{
					*cross_log << "\t" <<  costs[index_replace1] <<
						"\t" << costs[index_replace2] << std::endl;
				}
			/*	
            // mutate offspring, log mutation, update fitness
            // to be removed?

            // individual one
            double costTemp = costs[index_replace1]; 
            _organisms[index_replace1].mutate();
            if(costTemp != _organisms[index_replace1].getcost() && mut_log)
				{   
               *mut_log << _generation << "zzzz\t" << costs[index_replace1];
               costs[index_replace1] = costTemp;    
               *mut_log << "\t" << costs[index_replace1];
               costs[index_replace1] = _organisms[index_replace1].getcost(true);
               *mut_log << "\t" << costs[index_replace1] << std::endl;
            }
            // individual two
            costTemp = costs[index_replace2]; 
            _organisms[index_replace2].mutate();
            if(costTemp != _organisms[index_replace2].getcost() && mut_log)
				{   
               *mut_log << _generation << "zzzz\t" << costs[index_replace2];
               costs[index_replace2] = costTemp;    
               *mut_log << "\t" << costs[index_replace2];
               costs[index_replace2] = _organisms[index_replace1].getcost(true);
               *mut_log << "\t" << costs[index_replace2] << std::endl;
            }
            */
			}

			bool reeval = !(_generation % (costs.size()*10)) && _generation != 0;
			if(!(_generation % 1000) && _generation != 0)
			{
				// change light input source
				organism::swaplight();
				reeval = true;
			}
		
			// Force re-evaluation of all every now and then
			if(reeval)
			{
				for(size_t i = 0; i< costs.size(); ++i)
					costs[i] = _organisms[i].getcost(true);
			}

			if(log && log_every > 0 && (!(_generation % log_every) ||
				_generation + 1 == endgeneration))
			{
				int index_best = std::distance(costs.begin(),
					std::min_element(costs.begin(), costs.end()));

//test
/*
{
std::vector<size_t> numbers = _organisms[0].get_diversity();
for(size_t t : numbers)
   std::cout  << t << "\t";
std::cout << "\n";
}*/
				log_line(*log, _organisms[index_best], false);
			}

			if(savelog && save_every > 0 && !(_generation % save_every))
			{
				// save current state of best organism
				int index_best = std::distance(costs.begin(),
					std::min_element(costs.begin(), costs.end()));
				savelog->precision(20);
				*savelog << "generation: " << _generation <<  "\n";
				_organisms[index_best].save(*savelog);
			}

			_cputime += getcputime();
			if(_cputime > next_save)
			{
				save();
				next_save = _cputime + time_between_saves;
			}
		}
		
		save();
	}

	void anneal(size_t iterations, double t_start, double t_stop,
					std::ostream *log, unsigned save_every = 0, unsigned log_every = 0)
	{
		simulator::log = nullptr;
		assert(_organisms.size() >= 2);
		getcputime();

		if(log)
			log_header_anneal(*log);

		int index_cur = 0;  // current state
		int index_sug = 1;  // suggested state

		for(size_t i = 0; i < iterations; ++i, ++_generation)
		{

			double T = t_start * pow(t_stop / t_start, double(i) / iterations);

			// replace worst with our suggested update
			_organisms[index_sug] = _organisms[index_cur];
			_organisms[index_sug].mutate();

			double cost_cur = _organisms[index_cur].getcost();
			double cost_sug = _organisms[index_sug].getcost();

			// minimize cost
			bool accept = cost_sug <= cost_cur ||
				exp((cost_cur - cost_sug) / T) > gsl_rng_uniform(rng);

			if(accept)
				std::swap(index_cur, index_sug);

			if(log && log_every > 0 && (!(i % log_every) || i + 1 == iterations))
			{
				log_line(*log, T, _organisms[index_cur], false);
			}
			if(save_every > 0 && !(i % save_every))
				save();
		}
		_cputime += getcputime();
	}


	void debug_print(const std::vector<size_t> &m_in, const std::vector<size_t> &m_out)
	{
		//Debug info:
		std::cout << "#"<< m_in.size() << "\t" << m_out.size() << std::endl;
		if(m_in.size() > m_out.size())
		{
			for(size_t i = 0; i < m_out.size(); i++)
				std::cout << m_in[i] << "\t" << m_out[i] << std::endl;
			for(size_t i = m_out.size(); i < m_in.size(); i++)
				std::cout << m_in[i] << std::endl;
		}
		else
		{
			for(size_t i = 0; i < m_in.size(); i++)
				std::cout << m_in[i] << "\t" << m_out[i] << std::endl;
			for(size_t i = m_in.size(); i < m_out.size(); i++)
				std::cout << "\t" << "\t" << m_out[i] << std::endl;
		}

	}

	// evolve a network with structure specified by: m_in, m_out.
	void anneal_null(const organism &like, int iterations,
						  double t_start, double t_stop, std::ostream *log,
						  unsigned save_every = 0, unsigned log_every = 0)
	{
		// get that organism's characteristics in form of number of
		// inputs and outputs on each node
		getcputime();
		std::vector<size_t> m_in, m_out;
		like.get_degree_distribution(m_in, m_out);

		assert(_organisms.size() >= 2);

		if(log)
			log_header_anneal(*log);

		int index_cur = 0;  // current state
		int index_sug = 1;  // suggested state
		int indx = 0;       // count how many failed annealings
		std::vector<double> costs(2);
		do{
			_organisms[index_cur] = organism(organism::get_genome_length_init());
			for(size_t i = 0; i < costs.size(); ++i)
				costs[i]=_organisms[i].getcost_null(m_in, m_out);

			for(size_t i = _generation; i < iterations+_generation; ++i)
			{
				double T = t_start * pow(t_stop / t_start, double(i) / iterations);
				double cost_cur = costs[index_cur];

				// replace worst with our suggested update
				_organisms[index_sug] = _organisms[index_cur];
				_organisms[index_sug].mutate();

				double cost_sug = _organisms[index_sug].getcost_null(m_in, m_out);
				costs[index_sug] = cost_sug;

				// minimize cost
				bool accept = cost_sug <= cost_cur ||
					exp((cost_cur - cost_sug) / T) > gsl_rng_uniform(rng);

				if(accept){
					_organisms[index_cur] = _organisms[index_sug];
					costs[index_cur] = cost_sug;
				}
				if(log && log_every > 0 &&
					(!(i % log_every) || i + 1 == iterations + _generation || costs[index_cur] == 0))
				{
					log_line(*log, T, _organisms[index_cur], false);
				}
				if(save_every > 0 && !(i % save_every)) // save even if not found fitness=0 yet
					save();
				if(cost_sug == 0)
					break;
			}
/*			if(indx)
			{
				std::cout << "try:\t" << indx << "\t" << costs[index_cur] << std::endl;
				std::vector<size_t> m2_out, m2_in;
				_organisms[index_cur].get_degree_distribution(m2_in, m2_out);
				debug_print(m_in, m2_in);
				debug_print(m_out, m2_out);
			}*/
			indx++;
		}while(costs[index_cur] != 0);

		// Discard all but the best
		if(index_cur > 0)
			_organisms[0] = _organisms[index_cur];
		resize(1);
		
	
		_generation += iterations;
		_cputime += getcputime();
	}


	void setdirectory(const std::string &newdir)
	{
		_directory = newdir;
	}

	void save() const
	{
		mkdir(_directory.c_str(), 0777);
		std::string savefile = _directory + "/out.save";
		std::string savetmp = savefile + ".tmp";
		std::ofstream out(savetmp);
		if(!out)
			throw std::string("Error: unable to open save file " + savetmp);

		out.precision(20);

		out << "generation: " << _generation << "\n"
			"cputime: " << _cputime << "\n";

		out << "organisms-size: " << _organisms.size() << "\n";
		for(size_t i = 0; i < _organisms.size(); ++i){
			out << "organism: " << i << "\n";
			_organisms[i].save(out);
		}

		if(std::rename(savetmp.c_str(), savefile.c_str()))
			throw "failed to rename " + savetmp + " to " + savefile;
	}

	// Returns false if no such file
	bool load()
	{
		std::string loadpath = _directory + "/out.save";

		std::ifstream input(loadpath);
		if(!input)
			return false;

		std::string s;
		if(!(input >> s) || s != "generation:")
			throw std::string("Error when loading generation.") + s;
		input >> _generation;

		if(!(input >> s) || s != "cputime:")
			throw std::string("Error when loading cputime.") + s;
		input >> _cputime;

		if(!(input >> s) || s != "organisms-size:")
			throw std::string("Error when loading organisms-size.");
		size_t n_organisms;
		input >> n_organisms;

		_organisms.clear();
		for(size_t i = 0; i < n_organisms; ++i)
		{
			size_t j;
			if((input >> s) && s == "organism:")
			{
				if((input >> j) && j == i)
					_organisms.push_back(organism(input));
				else
					throw std::string("Error when loading an individual organism.");
			}
			else
				throw std::string("Error when loading organism.");
		}

		return true;
	}

	// Returns false if no such file
	bool loadbest(std::ostream &log, double prune_limit = -1)
	{
		std::string loadpath = _directory + "/best.save";

		std::ifstream input(loadpath);
		if(!input)
			return false;

		log_header(log, true);

		_organisms.clear();
		std::string s;
		while((input >> s) && s == "generation:")
		{
			input >> _generation;
			organism best(input);

			double cost;
			if(prune_limit >= 0)
				cost = best.prune(prune_limit);
			else
				cost = best.getcost(true);
			log_line(log, cost, best, true);
		}
		return true;
	}

	organism &findbest(double &best_cost)
	{
		size_t best_index = 0;
		best_cost = 1e6;

		for(size_t i = 0; i < _organisms.size(); ++i)
		{
			double cost = _organisms[i].getcost();
			if(best_cost > cost)
			{
				best_cost = cost;
				best_index = i;
			}
		}
		return _organisms[best_index];
	}

	void eval(std::ostream &log)
	{
		double cost;
		organism &best = findbest(cost);
//		best.prune(.005);
//		best.prune(.005);
		cost = best.getcost(true);

		log_header(log, true);
		log_line(log, best, true);
	}

private:

	// choose random individual
	int choose_tournament(const std::vector<double> &costs)
	{
		assert(costs.size() > 1);
		int index1 = gsl_rng_uniform_int(rng, costs.size());
		int index2 = gsl_rng_uniform_int(rng, costs.size() - 1);
		if(index2 >= index1)
			++index2;
		if(costs[index1] < costs[index2])
			return index1;
		return index2;
	}
   
   // returns the population diversity based on the hamming distance
   // between gene pairs in between organisms
   double get_pop_diversity()
   {
      double nonorm = 0;
      double totgenes = 0; 
      for(size_t i = 0; i < _organisms.size(); i++)
      {
         const std::vector<gene> &genes = _organisms[i].get_genes();

         for(size_t j = 0; j < genes.size(); j++)   
         {
            const gene &headgene = genes[j];
            totgenes += _organisms.size() - 1; 
            for(size_t k = 0; k < _organisms.size(); k++)
            {
               // don't compare within host organism 
               if(k == i)
                  continue;
               const std::vector<gene> &compgenes = _organisms[k].get_genes(); 
               double min_dist = std::numeric_limits<double>::max();
               
               for(size_t l = 0; l < compgenes.size(); l++)
               {
                  int dist = headgene.distance(compgenes[l]);
                  if(dist < min_dist)
                     min_dist = dist; 
               }
               nonorm += min_dist;
            }
         }
      }
      return nonorm/totgenes/gene::gene_size;
   }
   
   std::string _directory;
	std::vector<organism> _organisms;
	size_t _generation;                 // save iterations done
	double _cputime;							// save accumulated cpu time spent
};


#endif
