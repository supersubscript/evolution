#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <functional>

#include "genes.h"
#include "organism.h"
#include "population.h"
#include "boolean.h"
#include "binom.h"

gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);


void get_dist_of_maxinputs(std::string loadpath, std::map<int,int> &maxinputdist)
{
	std::string source_file = "unprunedrules.0.out";
	std::string f = loadpath + "/" + source_file;
	std::ifstream infile(f.c_str());

	std::string a;
	getline(infile, a); // throw away a line x 2
	getline(infile, a);

	if(infile)
	{
		while(getline(infile, a))
		{
			int inputs = 0, freq = 0;
			std::istringstream ss(a);
			ss >> inputs;
			ss >> freq;
			maxinputdist[inputs] = freq;
		}
	}
	else
		throw std::string("must run stat load to generate "+ source_file +" file first");
}


class rulecount
{
public:
	static const int small_rule_limit = 2;

	rulecount()
	{
		for(int i = 0; i <= small_rule_limit; i++)
		{
			// Make room for counting 2**2**i rules with i inputs
			_small_rule_count.push_back(
				std::vector<size_t>(1 << (1 << i), 0));
		}
	}

	// for very large rules that we can't build a truth table for
	void add_unpruned_rule(int inputsize)
	{
		_degree_distr_unpruned[inputsize]++;
	}

	//if small rule
	void add_rule(size_t inputsize, unsigned input)
	{
		assert(inputsize < _small_rule_count.size());
		assert(input < _small_rule_count[inputsize].size());

		// Merge equivalent rules such as A&!B and B&!A
		if(inputsize == 2 && (input & 6) == 4)
			input -= 2;

		_small_rule_count[inputsize][input]++;
		_degree_distr[inputsize]++;

		if(inputsize >= 1 && inputsize <= 2)
		{
			bool noncanal = inputsize == 2 && (input == 6 || input == 10);
			_canal_count[inputsize][noncanal ? 0 : inputsize]++;
		}
	}

	//if not small rule, use _canal_count
	void add_rule_canal(int inputsize, int inp_canal)
	{
		_canal_count[inputsize][inp_canal]++;
		_degree_distr[inputsize]++;
	}

	static const std::string &rulename(size_t inputs, size_t rule)
	{
		static const std::vector<std::vector<std::string>> names{
			{"FALSE", "TRUE"},                      // 0 inputs
			{"FALSE", "!A", "A", "TRUE"},           // 1 input
			{"FALSE", "!A & !B", "A & !B", "!A",    // 2 inputs
			"blaj", "blaj", "A != B", "!A | !B",
			"A & B", "A = B", "A", "A | !B",
			"blaj", "blaj", "A | B", "TRUE" }};
		return names[inputs][rule];
	}

	void print_unpruned_rules(std::string outfile) const
	{
		std::ofstream out(outfile.c_str());
		out << "#degree distribution before pruning\n";
		out << "#inputs\tcount\tfreq\n";
		size_t grand_total = 0;
		for(auto& deg : _degree_distr_unpruned)
			grand_total += deg.second;
		for(auto& deg : _degree_distr_unpruned)
		{
			out << deg.first << "\t" << deg.second << "\t" <<
				deg.second / double(grand_total) << "\n";
		}
		out.close();
	}

	void print_rules(std::ostream &out) const
	{
		out << "#degree distribution after pruning\n";
		out << "#inputs\tcount\tfreq\n";
		size_t grand_total = 0;
		for(auto& deg : _degree_distr)
			grand_total += deg.second;
		for(auto& deg : _degree_distr)
		{
			out << deg.first << "\t" << deg.second << "\t" <<
				deg.second / double(grand_total) << "\n";
		}
		out << "\n";

		for(int inputs = 0; inputs <= small_rule_limit; inputs++)
		{
			// Skip input sizes with no data
			auto it = _degree_distr.find(inputs);
			if(it == _degree_distr.end())
				continue;
			double total = it->second;

			auto &counts = _small_rule_count[inputs];
			out << "# inputs: " << inputs << "\n";
			out << "#rule\tcount\tfreq\tp_lower(P<0.01)\tp_upper(P<0.01)\n";

			for(size_t i = 0; i < counts.size(); i++)
			{
//				// Only print rules that actually appear
//				if(!counts[i])
//					continue;
				if(inputs == 2 && (i & 6) == 4)
					continue;

				auto lim = binom_p_limits(counts[i], total, .005);

				out << rulename(inputs, i) << "\t" << counts[i] <<
					"\t" << counts[i] / total << "\t"
					 << lim.first << "\t" << lim.second << "\n";
			}
			out << "\n";
		}

		for(auto& inputs : _canal_count)
		{
			// Skip input counts with no rules
			if(inputs.second.empty())
				continue;

			// We can't use [] because of const-ness but we know
			// inputs is a key in _degree_distr
			double total = _degree_distr.find(inputs.first)->second;

			out << "# inputs: " << inputs.first << "\n";
			out << "# canalizing\tcount\tfreq\n";

			for(auto& inp : inputs.second)
			{
				out << inp.first << "\t"
					 << inp.second << "\t"
					 << inp.second / total << "\n";
			}
			out << std::endl;
		}

		out << "# fraction canalizing on >= 1 input\n";
		out << "# inputs\tfrac\tsem\n";
		for(auto inputs : _canal_count)
		{
			int tot = _degree_distr.find(inputs.first)->second;
			int c = tot - inputs.second[0];
			double frac = (double)c / tot;
			double sd = std::sqrt(c * (1. - frac)) / (tot - 1);
			out << inputs.first << "\t" << frac << "\t" << sd << "\n";
		}

		out << "\n# fraction nested canalizing (canalizing on all inputs)\n";
		out << "# inputs\tfrac\tsem\n";
		for(auto inputs : _canal_count)
		{
			int tot = _degree_distr.find(inputs.first)->second;
			int c = inputs.second[inputs.first];
			double frac = (double)c / tot;
			double sd = std::sqrt(c * (1. - frac)) / (tot - 1);
			out << inputs.first << "\t" << frac << "\t" << sd << "\n";
		}

		out << "\n# number of canalizing inputs\n";
		out << "# inputs\tmean\tsem\n";
		for(auto inputs : _canal_count)
		{
			double sum = 0, sum2 = 0, n = 0;
			for(auto& inp : inputs.second)
			{
				double v = inp.first, w = inp.second;
				n += w;
				sum += w * v;
				sum2 += w * v * v;
			}

			double mean = sum / n;
			double sd = std::sqrt(std::max(0., sum2 / n - sqr(mean)) / (n - 1));
			out << inputs.first << "\t" << mean << "\t" << sd << "\n";
		}

	}

	size_t get_degree_distr(int n)
	{
		return _degree_distr[n];
	}

private:
	std::vector<std::vector<size_t>> _small_rule_count;

	// inp_used -> ( inp_canal -> count )
	std::map<int, std::map<int, size_t>> _canal_count;

	// In-degree distribution; inputs -> number of such rules
	std::map<int, size_t> _degree_distr;
	std::map<int, size_t> _degree_distr_unpruned;
};


unsigned compactify_rule(const std::vector<int> &rule)
{
	unsigned compactrule = 0;
	for(unsigned i = 0; i < rule.size(); ++i)
		if(rule[i])
			compactrule |= 1 << i;
	return compactrule;
}

// rnapthreshold positive: actual threshold. negative: relative highest
// observed value among the 2**n (e.g. -.5 for half max).
void examinerules(const gene &cisreg, std::vector<double> tfexpression,
	rulecount &storerule, std::ostream *out, double rnapthreshold)
{
	std::set<size_t> tf_index;
	cisreg.getbindingtfs(tf_index);
	int bindings = tf_index.size();

	storerule.add_unpruned_rule(bindings);

	// upper cap, time waits for no man
	if(bindings > 20)
		return;

	std::vector<double> RNAPoutput;
	organism::computerule(cisreg, tfexpression, RNAPoutput);

	if(rnapthreshold < 0)
	{
		double t = 0;
		for(size_t j = 0; j < RNAPoutput.size(); ++j)
			t = std::max(RNAPoutput[j], t);
		rnapthreshold = -rnapthreshold * t;
	}

	if(out)
	{
//		*out << "#binding sites:\n";
//		cisreg.printbindingsites(*out);
//		cisreg.printrnapactivities(*out);

		*out << "\n#output levels:\n";
		for(auto v : RNAPoutput)
			*out << v << "\t" << v / rnapthreshold << "\n";
	}

	// Binarize the cisregion output to make a Boolean rule
	std::vector<int> unpruned_rule(RNAPoutput.size());
	for(size_t j = 0; j < RNAPoutput.size(); ++j)
		unpruned_rule[j] = RNAPoutput[j] > rnapthreshold;

	// remove unused rules
	std::vector<int> rule;
	int inp_unused;
	prune_unused_input(unpruned_rule, rule, inp_unused);
	int inp_used = bindings - inp_unused;

	assert(1U << inp_used == rule.size());

	if(inp_used > rulecount::small_rule_limit)
	{
		// Find the number of nested canalizing inputs, and the
		// rule for the remaning inputs.
		int inp_canal;
		std::vector<std::pair<bool, bool>> canal_inout;
		std::vector<int> other_rule;

		// get nested canalizing part
		check_canalization(rule, other_rule, inp_canal, canal_inout);
		int inp_other = inp_used - inp_canal;

		storerule.add_rule_canal(inp_used, inp_canal);

		if(out && false)
		{
			*out << "\n";
			if(inp_unused)
				*out << "unused inputs: " << inp_unused << "\n";
			if(!canal_inout.empty())
			{
				*out << "canalizing rule:";
				for(auto v : canal_inout)
					*out << "  " << (int)v.first << ":" << (int)v.second;
				if(other_rule.size() == 1)
					*out << "  " << other_rule.front();
				*out << "\n";
			}
			if(other_rule.size() > 1 || canal_inout.empty())
			{
				*out << "other rule:";
				if(inp_other <= rulecount::small_rule_limit)
				{
					*out << "  " << rulecount::rulename(
						inp_other, compactify_rule(other_rule));
				}
				else
				{
					for(auto v : other_rule)
						*out << "  " << v;
				}
				*out << "\n";
			}
			*out << "\n";
		}
	}
	else
	{
		// Convert rule representation from vector or 0,1 to bit
		unsigned compactrule = compactify_rule(rule);
		storerule.add_rule(inp_used, compactrule);

		if(out && false)
			*out << "name: " << rulecount::rulename(inp_used, compactrule) << "\n\n";
	}
}


int printusage(char const* const* argv)
{
	std::string a = std::string(argv[0]);
	std::cerr <<
		"usage: " << a << " randomrules <inputs> <iterations> <tflevel> <rnaplevel> \n"
		"\t\t<promoter_start> <promoter_end> <enable_repressors>\n"
		"    or " << a << " randomrulesdist <loaddir> <iterations> <tflevel> <rnaplevel>\n"
		"    or " << a << " randombindingtype <loaddir> <iterations> \n"
		"    or " << a << " randomsamebinding <loaddir> <iterations>\n"
		"    or " << a << " savenet <loaddir> [dot (default) | tlp | sif | csv]\n"
		"    or " << a << " load <basepath> [optional tflevel]\n";
	return 1;
}


int main(int argc, char *argv[])
{
	gsl_rng_set(rng, randrandseed());

	if(argc < 2)
		return printusage(argv);

	// Extract/remember the directory we're running from
	sourcedir(argv[0]);

	std::string mode(argv[1]);
	int args;			// correct length of argc for
	int optargs = 0;  // allow some optional arguments

	if(mode == "randomrules")
		args = 9;
	else if(mode == "randomrulesdist")
		args = 6;
	else if(mode == "randomsamebinding")
		// for generating single bindings to random cis regions, and
		// check how many binding sites 0..n, and how many pos/neg.
		args = 4;
	else if(mode == "load")
	{
		// Process and gather all kind of fun statisitics on the
		// finished networks in the <basepath>
		args = 3;
		optargs = 1;
	}
	else if(mode == "savenet")
	{
		// Create an output file in some graph format for network in
		// <loaddir> provided.
		args = 3;
		optargs = 1;
	}
	else if(mode ==  "randombindingtype")
	{
		// Compute probability of TF-DNA interaction bindingtypes,
		// (competetive, cooperative, or other) expected from a
		// random/neutral net.
		args = 2+2;
		optargs = 0;
	}
	else
		return printusage(argv);

	if(argc < args || argc > args + optargs)
		return printusage(argv);

	try
	{
		if(mode == "randomrules")
		{
			int inputs = toint(argv[2]);
			size_t iterations = tosize_t(argv[3]);
			double tflevel = todouble(argv[4]);
			double rnapmaxoutput = todouble(argv[5]);

			int pstart = toint(argv[6]);
			int pend = toint(argv[7]);
			bool repressors = toint(argv[8]);

			gene::set_promoter(pstart, pend);
			gene::set_enable_repressors(repressors);

			rulecount rules; // stores rule statistics

			// create a cisregion with $input TFs
			//loop for statistics
			while(rules.get_degree_distr(2) < iterations)
			{
				gene cis;
				cis.randomize();
				cis.randomize_sites(inputs);
				std::vector<double> tflevels(inputs, tflevel);
				examinerules(cis, tflevels, rules, nullptr, rnapmaxoutput);
			}

			mkdir("random", 0777);
			std::ostringstream fname;
			fname << "random/rules." << pstart << "-" << pend << ".r" <<
				repressors << ".tf" << tflevel << ".rnap" << rnapmaxoutput <<
				".in" << inputs << ".it" << iterations << ".out";
			std::ofstream out(fname.str());

			out << "# iterations " << iterations << "\n"
				"# tflevel " << tflevel << "\n\n";

			rules.print_rules(out);
		}
		if(mode == "randomrulesdist")
		{
			std::string loadpath = std::string(argv[2]);
			int iterations = toint(argv[3]);
			double tflevel = todouble(argv[4]);
			double rnapmaxoutput = todouble(argv[5]);
			organism::load_globals(loadpath);

			rulecount rules; // stores rule statistics

			// read in nunber or maxinputs from a distribution from the
			// network we are to load
			std::map<int, int> maxinputdist;
			get_dist_of_maxinputs(loadpath, maxinputdist);

			int number_of_rules = 0;
			do{
				for(auto &inp : maxinputdist)
				{
					int inputs = inp.first;

					//loop for statistics
					for(int i = 0; i < inp.second; i++)
					{
						gene cis;
						cis.randomize();
						cis.randomize_sites(inputs);
						std::vector<double> tflevels(inputs, tflevel);
						examinerules(cis, tflevels, rules, nullptr, rnapmaxoutput);
					}
					number_of_rules += inp.second;
				}
			}while(number_of_rules < iterations);

			std::ostringstream fname;
			fname << loadpath << "/rules" << ".tf"
					<< tflevel << ".it" << iterations << ".out";
			std::ofstream out(fname.str());

			out << "# iterations " << iterations << "\n"
				"# tflevel " << tflevel << "\n\n";

			rules.print_rules(out);
		}
		if(mode == "randomsamebinding")
		{
			// for generating single bindings to random cis regions, and
			// check how many binding sites 0..n, and how many pos/neg.

			std::string loadpath = std::string(argv[2]);
			int iterations = toint(argv[3]);
			organism::load_globals(loadpath);
			std::map<std::pair<int, int>, int> hist_sign_same;

			for(int i = 0; i < iterations ; ++i)
			{
				gene cis;
				cis.randomize();
				std::map<int, size_t> stump_hist_binding_same;         // not used
				std::map<size_t, std::vector<std::pair<int, int>>> hist_sign;

				cis.try_binding(1,stump_hist_binding_same);
				const size_t genes_in_organism = 1;

				cis.get_tfstat(stump_hist_binding_same, hist_sign, genes_in_organism);
				hist_sign_same[hist_sign[0][0]]++;
			}
			std::ostringstream fname;
			fname << loadpath << "/randsamebinding"
					<< ".it" << iterations << ".out";
			std::ofstream out(fname.str());

			out << "# iterations " << iterations << "\n"
				 << "# sign of TFs binding to same gene\n" << "pos\tneg\tfreq" << std::endl;
			for(auto& bindsign : hist_sign_same)
				out << bindsign.first.first << "\t" << bindsign.first.second
					 << "\t" << bindsign.second << std::endl;
				out.close();
		}
		else if(mode == "savenet")
		{
			std::string loadpath = std::string(argv[2]);
			std::string format = args + optargs == argc ? std::string(argv[3]) : "dot";

			if(format != "csv" && format != "tlp" && format != "sif" && format != "dot")
				throw "format must be \'cvs\', \'tlp\', \'dot\', or \'sif\' not: " + format;

			// Check if loadpath is a directory with many networks;
			// otherwise assume it's a single one.
			std::set<std::string> dirnames;
			std::ifstream infile(std::string(loadpath + "/records"));
			if(infile)
			{
				// find which networks (populations) are done
				std::string a;
				while(getline(infile, a))
				{
					std::istringstream ss(a);
					ss >> a;
					dirnames.insert(loadpath + "/" + a);
					organism::load_globals(loadpath);
				}
			}
			else{
				organism::load_globals(loadpath + "/..");
				dirnames.insert(loadpath);
			}

			for(auto dirname : dirnames)
			{
				// load organisms
				population pop(dirname);
				if(!pop.load())
					throw "Failed to load from " + dirname;

				std::string savefile = dirname + "/network" + "." + format;
				std::cout << "saving net to\t" << savefile << std::endl;

				// find best organism
				double best_cost;
				organism &best = pop.findbest(best_cost);
				best.savenetwork(savefile, format);
			}
		}
		else if(mode == "load")
		{
			double tflevel = argc > 3 ? todouble(argv[3]) : 0;
			std::string loadpath = std::string(argv[2]);

			rulecount rules; // stores rule statistics

			// Assume loadpath is a directory with many networks
			organism::load_globals(loadpath);
			std::set<std::string> dirnames;
			std::ifstream infile(std::string(loadpath + "/records"));
			if(!infile)
				throw std::string("failed to open records file");
			// find which networks (populations) are done
			std::string a;
			while(getline(infile, a))
			{
				std::istringstream ss(a);
				ss >> a;
				dirnames.insert(loadpath + "/" + a);
			}

			std::string savefit = loadpath + "/fitness.hist";
			std::ofstream fitnessout(savefit.c_str());
			fitnessout << "# fitness" << std::endl;

			std::map<int, size_t> hammingHist;    // histogram, hamming dist of bindingsites
			std::map<int, size_t> bindSameHist;   // histogram, Tf binding to same gene
			std::map<std::pair<int, int>, int> bindSameHistSign;  // histogram, sign of Tf binding to same gene
			std::map<std::pair<int, int>, int> bindGoodHistSign;  // ditto, for TFs with ~50% positive binding
			std::map<size_t, std::vector<size_t>> bindPosNegHist; // histogram, pos/neg freq of bindingsites

			// count pairwise characteristics of binding sites overlap.
			// Keep in vector to get standard deviation
			typedef organism::bindingstats bindingstats;
			std::vector<bindingstats> samebindingstats, otherbindingstats;

			// store where on cis region bindingsites start
			std::map<size_t, size_t> binding_pos;

			// save #nodes, edges, in_degree, out_degree for best org.
			// for each net. (for bug hunting, and code checking)
			std::string degree_path = loadpath + "/degree_dist";
			std::ofstream degree_distributions(degree_path);

			std::vector<int> cisonebits(gene::cisregion_size + 1);
			std::vector<int> tfonebits(gene::bindmotif_size + 1);

			std::ofstream expressionlevels_out;
			if(!tflevel)
			{
				std::string fname = loadpath + "/expressionlevels";
				expressionlevels_out.open(fname);
				if(!expressionlevels_out)
					throw "unable to open " + fname;
				expressionlevels_out << "#expressionlevel\trnapactivity\tdecayrate\n";
			}

			for(auto dirname : dirnames)
			{
				std::cout << "#doing:\t" << dirname << std::endl;

				// load organisms
				population pop(dirname);
				if(!pop.load())
					throw "Failed to load from " + dirname;

				// find best organism
				double best_cost;
				organism &best = pop.findbest(best_cost);
				best.printbehaviour(dirname);

				fitnessout << best_cost << std::endl;
				// get its cisregions
				std::vector<gene> cisregs;
				best.get_cisregions(cisregs);

				// The 'on' level of TFs, from dynamics or command line
				std::vector<double> maxexpression, maxrnapactivity, maxdecayrate;
				if(tflevel)
					maxexpression.resize(cisregs.size(), tflevel);
				else
				{
					best.findmaximumexpression(maxexpression, maxrnapactivity, maxdecayrate);
					for(size_t i = 0; i < maxexpression.size(); ++i)
					{
						expressionlevels_out << maxexpression[i] << "\t" <<
							maxrnapactivity[i] << "\t" << maxdecayrate[i] << "\n";
					}
					expressionlevels_out << std::endl;
				}

continue;

				// loop through its cisregs
				for(size_t i = 0; i < cisregs.size(); ++i)
				{
					examinerules(cisregs[i], maxexpression, rules, nullptr,
						tflevel ? -.5 : .5 * maxrnapactivity[i]);
				}

				// save info for histogram of hamming distances
				best.gethammingdistances(hammingHist);

				// Get some TF stats
				bindingstats samebinding, otherbinding;
				best.get_tfstats(bindSameHist, bindSameHistSign, bindGoodHistSign,
					bindPosNegHist, samebinding, otherbinding);
				samebindingstats.push_back(samebinding);
				otherbindingstats.push_back(otherbinding);

				// get some stat on in/out distribution
				degree_distributions << dirname << std::endl;
				best.get_degreestat(degree_distributions);

				// get stat on where on a cis the binding site is located
				best.get_bindingposstat(binding_pos);

				// get stat on number of 1-bits
				best.get_onebitstat(cisonebits, tfonebits);
			}
			fitnessout.close();
			degree_distributions.close();

			// Place the results together with the networks
			std::ostringstream fname1, fname2;
			fname1 << loadpath << "/rules." << tflevel << ".out";
			fname2 << loadpath << "/unprunedrules." << tflevel << ".out";
			std::cout << "#saving rule stats to " << fname1.str() << "\n";

			rules.print_unpruned_rules(fname2.str());

			std::ofstream out(fname1.str());
			if(tflevel)
				out << "# tflevel = " << tflevel << "\n\n";
			else
				out << "# tflevel = dynamic_max\n\n";
			rules.print_rules(out);
			out.close();

			// Make a histogram of hamming distances for all bindings,
			// ------------------------------------------------------
			{
				// for the best organism in each net
				std::string hammingpath = loadpath + "/hamming.hist";
				std::ofstream out_hamming(hammingpath.c_str());

				out_hamming << "# hamming\t count" << std::endl;
				for(auto &hamming : hammingHist)
					out_hamming << hamming.first << "\t" << hamming.second << std::endl;
				out_hamming.close();
			}
			// Make a histogram of tf binding to same gene
			// -------------------------------------------
			{
				std::string savefile_bind = loadpath + "/tf_binding.hist";
				std::cout << "#saving tf binding histogram to\t" << savefile_bind << std::endl;
				std::ofstream out_bind(savefile_bind);

				out_bind << "#tf_bind_to_same\t freq" << std::endl;
				for(auto& bind : bindSameHist)
					out_bind << bind.first << "\t" << bind.second << std::endl;
				out_bind.close();
			}

			// Make a histogram of sign of tf binding to same gene (sign ambiguity)
			// --------------------------------------------------------------------
			{
				std::string savefile_bind_sign = loadpath + "/tf_posneg_same.hist";
				std::cout << "#saving same tf sign binding histogram to\t" << savefile_bind_sign << std::endl;
				std::ofstream out_bind_sign(savefile_bind_sign);

				out_bind_sign << "#sign of TFs binding to same gene\n" << "pos\tneg\tfreq" << std::endl;
				for(auto& bindsign : bindSameHistSign)
					out_bind_sign << bindsign.first.first << "\t" << bindsign.first.second
									  << "\t" << bindsign.second << std::endl;
				out_bind_sign.close();
			}

			// Make a histogram of sign of tf binding to same gene (sign ambiguity)
			// --------------------------------------------------------------------
			{
				std::string savefile_bind_sign = loadpath + "/tf_posneg_good.hist";
				std::cout << "#saving same tf sign binding histogram to\t" << savefile_bind_sign << std::endl;
				std::ofstream out_bind_sign(savefile_bind_sign);

				out_bind_sign << "#sign of TFs binding to same gene\n" << "pos\tneg\tfreq" << std::endl;
				for(auto& bindsign : bindGoodHistSign)
					out_bind_sign << bindsign.first.first << "\t" << bindsign.first.second
									  << "\t" << bindsign.second << std::endl;
				out_bind_sign.close();
			}

			// Make file with data on binding types
			// ------------------------------------------------------
			{
				auto f = [](const std::vector<bindingstats> &values, double bindingstats::*member,
					bindingstats &meanb, bindingstats &stderrb)
				{
					double mean = 0, stderr = 0;
					double N = values.size();

					for(auto& value : values)          // compute mean value
						mean += value.*member;
					mean /= N;
					meanb.*member = mean;

					for(auto& value : values)          // compute standard dev.
						stderr += sqr(mean - value.*member);
					stderr /= N * (N-1);
					stderrb.*member = sqrt(stderr);
				};

				auto g = [&f, &loadpath](const std::vector<bindingstats> &stats, const std::string &type)
				{
					std::string bindtype = loadpath + "/bindtype_" + type;
					std::ofstream out_bindtype(bindtype.c_str());

					bindingstats mean, stderr;

					f(stats, &bindingstats::overlapping, mean, stderr);
					f(stats, &bindingstats::cooperative, mean, stderr);
					f(stats, &bindingstats::other, mean, stderr);

					double sum = mean.overlapping + mean.cooperative + mean.other;

					out_bindtype
						<< "#TYPE\t MEAN\t ERR\n"
						<< "Overlap\t" << mean.overlapping / sum
						<< "\t" << stderr.overlapping / sum
						<< "\nCooperative\t" << mean.cooperative / sum
						<< "\t" << stderr.cooperative / sum
						<< "\nOther\t" << mean.other / sum
						<< "\t" << stderr.other / sum
						<< std::endl;
				};

				g(samebindingstats, "same");
				g(otherbindingstats, "other");
			}
			// Make a histogram of tf pos/neg bindingsite
			// ------------------------------------------------------
			{
				std::string savefile_posneg = loadpath + "/tf_posneg.hist";
				std::cout << "#saving tf posneg histogram to\t" << savefile_posneg << std::endl;
				std::ofstream out_posneg(savefile_posneg);

				out_posneg << "#pos+neg\t pos_bindings\t freq" << std::endl;
				for(auto& sum : bindPosNegHist)
				{
					for(size_t i = 0; i <  sum.second.size(); i++)
						out_posneg << sum.first << "\t" << i << "\t" << sum.second[i] << "\n";
					out_posneg << std::endl;
				}
				out_posneg.close();
			}
			// Output data on where bindingsites start on cis regions
			// ------------------------------------------------------
			{
				std::string savefile_bindpos = loadpath + "/bindpos";
				std::cout << "#saving binding position to\t" << savefile_bindpos << std::endl;
				std::ofstream out(savefile_bindpos);
				out << "#site\tfreq" << std::endl;
				for(auto& pos : binding_pos)
					out << pos.first << "\t" << pos.second << std::endl;
				out.close();
			}

			{
				std::string savefile = loadpath + "/cisonebits.hist";
				std::ofstream out(savefile);
				out << "#onebits\t freq" << std::endl;
				for(size_t i = 0; i < cisonebits.size(); ++i)
					out << i << "\t" << cisonebits[i] << std::endl;
			}
			{
				std::string savefile = loadpath + "/tfonebits.hist";
				std::ofstream out(savefile);
				out << "#onebits\t freq" << std::endl;
				for(size_t i = 0; i < tfonebits.size(); ++i)
					out << i << "\t" << tfonebits[i] << std::endl;
			}

		}
		else if(mode == "randombindingtype")
		{
			// Read in settings file of path given

			std::string loadpath = std::string(argv[2]);
			int iterations = toint(argv[3]);
			organism::load_globals(loadpath);

			typedef organism::bindingstats bindingstats;
//			std::vector<bindingstats> samebindingstats, otherbindingstats;

			std::map<int, size_t> bindSameHist;   // histogram, Tf binding to same gene
			std::vector<int> cisonebits(gene::cisregion_size + 1);
//			std::vector<int> tfonebits(tfmotif_t::bits + 1);

			bindingstats ssum, ssum2, osum, osum2;

			auto ff = [](const bindingstats &value, bindingstats &sum, bindingstats &sum2)
			{
				sum.overlapping += value.overlapping;
				sum.cooperative += value.cooperative;
				sum.other += value.other;
				sum2.overlapping += sqr(value.overlapping);
				sum2.cooperative += sqr(value.cooperative);
				sum2.other += sqr(value.other);
			};

			for(int i = 0; i < iterations; i++)
			{
				gene cis;
				cis.randomize();

				// Try to create one pair of TFs for the cis region (might
				// come out with 0,1,2 binding sites)
				cis.try_binding(2, bindSameHist);

				// Check binding type of each pair of sites
				bindingstats samebinding, otherbinding;
				cis.get_bindingstat(samebinding, true);
				cis.get_bindingstat(otherbinding, false);

				ff(samebinding, ssum, ssum2);
				ff(otherbinding, osum, osum2);
//				samebindingstats.push_back(samebinding);
//				otherbindingstats.push_back(otherbinding);

				cisonebits[cis.onebits_cis()]++;
			}

			auto f = [](const bindingstats &sum, const bindingstats &sum2,
							double n, double bindingstats::*member,
							bindingstats &meanb, bindingstats &stderrb)
			{
				meanb.*member = sum.*member / n;
				double var = sum2.*member / n - sqr(meanb.*member);
				stderrb.*member = sqrt(std::max(0., var) / (n-1));
			};

			auto g = [&f, &loadpath, &iterations](const bindingstats &sum,
				const bindingstats &sum2, const std::string &type)
			{
				std::string bindtype = loadpath + "/bindtype_" + type;
				std::ofstream out_bindtype(bindtype.c_str());

				bindingstats mean, stderr;

				f(sum, sum2, iterations, &bindingstats::overlapping, mean, stderr);
				f(sum, sum2, iterations, &bindingstats::cooperative, mean, stderr);
				f(sum, sum2, iterations, &bindingstats::other, mean, stderr);

				double s = mean.overlapping + mean.cooperative + mean.other;

				out_bindtype
				<< "#TYPE\t MEAN\t ERR\n"
				<< "Overlap\t" << mean.overlapping / s
				<< "\t" << stderr.overlapping / s
				<< "\nCooperative\t" << mean.cooperative / s
				<< "\t" << stderr.cooperative / s
				<< "\nOther\t" << mean.other / s
				<< "\t" << stderr.other / s
				<< std::endl;
			};

			g(ssum, ssum2, "same_predicted");
			g(osum, osum2, "other_predicted");

			std::string savefile_bind = loadpath + "/tf_binding_predicted.hist";
			std::ofstream out_bind(savefile_bind);
			out_bind << "#tf_bind_to_same\t freq" << std::endl;
			for(auto& bind : bindSameHist)
				out_bind << bind.first << "\t" << bind.second << std::endl;

			std::string savefile_onebits = loadpath + "/cisonebits_predicted.hist";
			std::ofstream out_onebits(savefile_onebits);
			out_onebits << "#onebits\t freq" << std::endl;
			for(size_t i = 0; i < cisonebits.size(); ++i)
				out_onebits << i << "\t" << cisonebits[i] << std::endl;

		}
	}
	catch(std::string s)
	{
		std::cerr << "Error: " << s << "\n";
		return 1;
	}

	return 0;
}
