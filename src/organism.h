#ifndef DIGORG_ORGANISM_H
#define DIGORG_ORGANISM_H

#include <ostream>
#include <vector>
#include <memory>
#include <set>
#include "genes.h"
#include "bitstring.h"
#include "simulator.h"
#include "lightlevelgetter.h"

struct ode_parameters;



class organism
{
	// Cost function penalty on gene count
	static double _cost_penalty_genes;
	// Cost penalty on sum(squared in-degree)
	static double _cost_penalty_links;

	// Probability of insertion/deletion event
	static double _gene_indel_prob;
	// Power law controlling insertion/deletion length
	static double _indel_length_power;
	// Probability to flip a bit
	static double _gene_mutation_prob;
	// Probability to perform a crossover
	static double _gene_crossover_prob;
	// Probability of crossover length
	static double _gene_crossover_length_prob;

	// Max Hamming distance for matching 64-bit block in alignment divide-and-conquer
	static int _alignment_optimization_maxdiff;

	// Max length of our genome
	static size_t _genome_length_max;
	// Start length of our genome
	static size_t _genome_length_init;

public:
	enum cost_type
	{
		cost_undefined, cost_clock_ld, cost_clock_hf
	};

	static inline cost_type get_cost_type() { return _cost_type; }
	static void set_cost_type(cost_type v);

	// Loads from file
	explicit organism(std::istream &in);

	organism(size_t initial_genome_length);

	void reinitialize(size_t initial_genome_length);

	inline size_t size() const
	{
		return _genes.size();
	}

	inline size_t genome_size() const
	{
		return _genome.bits();
	}

	//mutate random gene(s), returning info about mutation types
	void mutate();

   const std::vector<gene>& get_genes() const
   {
      return _genes;
   }
     
	static inline double get_gene_mutation_prob()
	{
		return _gene_mutation_prob;
	}


	static inline double get_gene_crossover_prob()
	{
		return _gene_crossover_prob;
	}

	//genetic crossover
	static void crossover(const organism &source1, const organism &source2,
		organism &dest1, organism &dest2);

private:
	// Solve the ODE system, using GSL; store integral of protein expression
	// in expression_windows. The simulation runs from time 0 to time
	// expression_windows[d].size() * window_size.
	// expression_windows[d] has one element per time window. Each element
	// is a vector which will be resized by solve so that the levels of
	// all the proteins can be stored there.
	bool solve(std::vector<std::vector<std::vector<double>>> &expression_windows,
		size_t window_count, const lightlevelgetter &lightgetter,
		bool findmaxexpression, bool updateproteinlevels);

	// Entrains the clock system until some measure of convergence is
	// fulfilled or a number of days have passed.
	bool entrain(const lightlevelgetter &lightgetter);

	// For cost_clock, or rule for cost_majority
	double getcost_clock_ld(std::vector<size_t> &window_genes, bool fixed_output_genes);
	double getcost_clock_hf(std::vector<size_t> &window_genes, bool fixed_output_genes);

	// Computes the cost over a set of light conditions, after
	// entraining in those same conditions.
	double getcost(std::vector<size_t> &output_genes, bool fixed_output_genes);

public:
	// This one doesn't load/save info about output genes
	double getcost(bool forceEval = false);

	// costfuntion for a network to evolve into one with
	// architecture/structure defined by inputs
	double getcost_null(const std::vector<size_t> &m_in,
		const std::vector<size_t> &m_out) const;

	// Changes the light conditions used in the fitting
	static void swaplight();

	// Finds nodes that in one way or other affect the output node(s)
	void find_used_nodes(const std::vector<size_t> &output_genes,
		std::set<size_t> &used_genes) const;

	// Iteratively removes unused links until the cost has increased by no
	// more than the given amount. Returns new cost.
	double prune(double cost_increase_limit);

	// Returns the number of links in the network (TFs binding
	// cisregions, counting multiple binding sites in the same
	// cisregion only once, or the same thing for protein regulating
	// degraders/complex depending on the bool argument).
	int getlinkcount(bool regulators, bool degraders, const std::set<size_t> *used_genes = nullptr) const;

	// Returns the number that goes into the cost term for links. This is the
	// sum over cisregions of the cube of the number of TFs that bind that
	// cisregion, beyond the first two. That is, we only penalize cisregions
	// with three or more inputs.
	int getexcessivelinkcount() const;

	// Finds the min/max protein levels and the maximum transcription rates
	// in the current light conditions.
	void findmaximumexpression(std::vector<double> &maxexpression,
		std::vector<double> &maxrnapactivity, std::vector<double> &maxdecayrate);

	static void computerule(const gene &cisreg, const std::vector<double> &maxexpression,
		std::vector<double> &RNAPoutput);

	//save/load state to/from stream.
	void save(std::ostream &out) const;
	void load(std::istream &in);

	// investigate distribution of hamming disctances
	void gethammingdistances(std::map<int, size_t> &hamminghistogram) const;

	void get_degreestat(std::ostream &out) const;

	// return start site for bindingsites for statistics
	void get_bindingposstat(std::map<size_t,size_t> &positions);

	// return start site for bindingsites for statistics
	void get_onebitstat(std::vector<int> &cisonebits, std::vector<int> &tfonebits) const
	{
		for(size_t i = 0; i < _genes.size(); ++i)
		{
			cisonebits[_genes[i].onebits_cis()]++;
			tfonebits[_genes[i].onebits_motif()]++;
		}
	}

	// needed in stat to do rule statisitcs on the evolved cis-regions
	void get_cisregions(std::vector<gene> &cisreg) const
	{
		cisreg = _genes;
	}

	// Fills up the vectors with one element per tf/cisregion with the number
	// of inputs and outputs of that node, based on whether each TFs binding
	// to each cisregion. The vectors are sorted from low to high.
	void get_degree_distribution(std::vector<size_t> &m_in,
		std::vector<size_t> &m_out) const;

	static void set_genome_length_init(size_t l)
	{
		_genome_length_init = l;
	}
	static size_t get_genome_length_init()
	{
			return _genome_length_init;
	}

	typedef gene::bindingstats bindingstats;

	// Return a measure of how much the most oscillating gene oscillates.
	void oscillation(std::vector<double> &best_expr, std::vector<size_t> &best_index);

	// return number of all (inc. overlapping) paths which loop back,
	// defined either by TF-gene interactions, protein-protein interactions,
	// or both.
	size_t countloops(bool transcriptional, bool posttranslational) const;

	// investigate frequency of TF binding sites characteristics.
	void get_tfstats(std::map<int, size_t> &countSameBinding,
		std::map<std::pair<int, int>, int> &countSameBindingSign,
		std::map<std::pair<int, int>, int> &countGoodBindingSign,
		std::map<size_t, std::vector<size_t>> &poshist,
		bindingstats &samebinding, bindingstats &otherbinding) const;

	// print network in a firendly format (CVS)
	void savenetwork(std::string savefile, std::string format) const;

	// Prints timecourse to a set of files named by the light conditions.
	void printtimecourse(const std::string &basefname);

	// Prints timecourses
	void printbehaviour(const std::string &dirname);

	static void load_globals(const std::string &dir);
	static void save_globals(const std::string &dir);

	static void log_header(std::ostream &log, bool full);
	void log_line(std::ostream &log, bool full);
	std::vector<int> get_genetypes(const std::set<size_t> *used_genes) const;

	void log_cross(std::ostream &log, bool full);
   std::vector<size_t> get_diversity();
private:
	// Update all binding sites and (re)initialize _odesolver for
	// the right number of variables.
	void refresh();

	// For things that are common to the public constructors
	organism();

	//We need a static function in order to call compute(), since it's a member
	//function, but since "func" is static it behaves almost like a global
	//function, which is what we need.
	static int odefunc(double t, const double y[], double f[], void *params);

	//ode-system, on GSL-friendly form:
	int compute(double light, const double y[], double f[]);

	//fill vector with all expressionlevels (in each time window) for
	//all genes and light conditions
	bool getexpressionlevels(size_t window_count,
		std::vector<std::vector<std::vector<std::vector<double>>>> &all_expression_windows,
		const std::vector<std::shared_ptr<lightlevelgetter>> &light_conditions, bool norm);

	// Remember fitness when we know it, if unknown/outdated: -1
	double _cost;

	// Light conditions to use in the clock cost function right now
	static std::vector<std::shared_ptr<lightlevelgetter>> _light_conditions;

	// All possible light conditions available to us
	static std::vector<std::shared_ptr<lightlevelgetter>> _all_light_conditions;

	// Know which index in _all_light_conditions is next in line for
	// being used by _light_conditions
	static size_t _next_light;

	// The entire (modelled) DNA of the organism:
	bitstring _genome;

	// Store the genes from _genome in a vector
	std::vector<gene> _genes;

	// The type of cost function used for all organisms
	static cost_type _cost_type;

	//store gsl setup parameters (err_rel, abs. etc).
	mutable simulator::control _sim_control;

	//our GSL ODE-solver, as ponter so all different models/equ.systems
	//can use the same solver.
	std::shared_ptr<simulator> _odesolver;
};

#endif
