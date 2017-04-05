#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "common.h"
#include "bitstring.h"
#include "alignment.h"


size_t default_bslen = 10000;
int default_iters = 100;
static int heuristic_maxdiff = 20;
static size_t synapsing_minlen = 3;
static size_t synapsing_alignlen = 20;
static int genome_maxlen_factor = 2;

static enum { fit_orig, fit_barrier4, fit_barrier8 } fitness_function = fit_orig;

gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);



void print(std::ostream &out, const std::vector<std::vector<double>> &vecs)
{
	for(auto& vec : vecs)
	{
		double sum = 0, sum2 = 0, n = vec.size();
		for(auto v : vec)
		{
			sum += v;
			sum2 += v * v;
		}
		double mean = sum / n;
		double sd = std::sqrt(std::max(0., sum2 / n - sqr(mean)) / (n - 1));
		out << "\t" << mean << "\t" << sd;
	}
}

void print(std::ostream &out, double x,
	const std::vector<std::vector<std::vector<double>>> &vecs)
{
	out << x;
	for(auto &v : vecs)
		print(out, v);
	out << std::endl;
}

void prepare_history(bitstring &bs, std::vector<int> &history)
{
	history.resize(bs.bits());
	bool prev = bs[0];
	int j = 0;
	for(size_t i = 0; i < bs.bits(); ++i)
	{
		if(bs[i] != prev)
		{
			prev = bs[i];
			++j;
		}
		history[i] = j;
	}
}

void indel_with_history(bitstring &bs, std::vector<int> &history, double indel_length_power)
{
	// duplicate/insert AND delete
	size_t len = ran_discrete_power(rng, indel_length_power, bs.bits());
	size_t from = gsl_rng_uniform_int(rng, bs.bits() - len + 1);
	size_t to = gsl_rng_uniform_int(rng, bs.bits() + 1);
	bs.insert(to, bs.substr(from, len));
	history.reserve(history.size() + len);
	history.insert(history.begin() + to,
		history.begin() + from, history.begin() + from + len);
	size_t start = gsl_rng_uniform_int(rng, bs.bits() - len);
	bs.remove(start, len);
	history.erase(history.begin() + start, history.begin() + start + len);
}

int score_alignment_history(const std::vector<char> &avec, const std::vector<char> &bvec,
	std::vector<int> &ahist, std::vector<int> &bhist)
{
	int score = 0;
	for(size_t apos = 0, bpos = 0, i = 0; i < avec.size(); ++i)
	{
		if(avec[i] != '-' && bvec[i] != '-' && ahist[apos] == bhist[bpos])
			++score;
		if(avec[i] != '-')
			++apos;
		if(bvec[i] != '-')
			++bpos;
	}
	return score;
}


static std::ostream &operator<<(std::ostream &out, const std::vector<int> &vec)
{
	for(auto v : vec)
		out << "\t" << v;
	return out;
}

static void test_indel_history()
{
	size_t bslen = 100;
	double mutrate = 0;
	double indelrate = 0.02;
	double indel_length_power = -1.5;

	bitstring bs1(bslen);
	bs1.rand();
	std::vector<int> history1;
	prepare_history(bs1, history1);
	bitstring bs2(bs1);
	std::vector<int> history2(history1);

	int indels = (int)(indelrate * bslen + .5);
	double mrate = mutrate / (indels + 1);
	for(int j = 0; j < indels; ++j)
	{
		bs2.mutate(mrate);
		indel_with_history(bs2, history2, indel_length_power);
	}
	bs2.mutate(mrate);

	std::vector<char> a1, a2;
	align_hirschberg(bs1, bs2, a1, a2);
	for(auto c : a1)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\n";
	for(auto c : a2)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\n";
	std::cout << history1 << "\n" << history2 << "\n";
	std::cout << "score " << score_alignment_history(a1, a2, history1, history2) << "\n";
}

// data is time / score / newscore => (hirschberg / heuristic / synapsing => (iteration => value))
void time_methods(int iters, size_t bslen, double mutrate, double indelrate,
	std::vector<std::vector<std::vector<double>>> &data, bool skip_slow = false)
{
	double indel_length_power = -1.5;

	data.clear();
	data.resize(3, std::vector<std::vector<double>>(4));

	for(int i = 0; i < iters; ++i)
	{
		bitstring bs1(bslen);
		bs1.rand();
		std::vector<int> history1;
		prepare_history(bs1, history1);
		bitstring bs2(bs1);
		std::vector<int> history2(history1);

		int indels = (int)(indelrate * bslen + .5);
		double mrate = mutrate / (indels + 1);
		for(int j = 0; j < indels; ++j)
		{
			bs2.mutate(mrate);
			indel_with_history(bs2, history2, indel_length_power);
		}
		bs2.mutate(mrate);

		for(int j = skip_slow ? 2 : 0; j < 4; ++j)
		{
			std::vector<char> a1, a2;
			getcputime();
			if(!j)
				align_hirschberg(bs1, bs2, a1, a2);
			else if(j == 1)
				align_heuristic(bs1, bs2, a1, a2, heuristic_maxdiff);
			else if(j == 2)
				align_synapsing(bs1, bs2, a1, a2, synapsing_minlen, 0);
			else
				align_synapsing(bs1, bs2, a1, a2, synapsing_minlen, synapsing_alignlen);

			data[0][j].push_back(getcputime());
			data[1][j].push_back(score_alignment(a1, a2));
			data[2][j].push_back(score_alignment_history(a1, a2, history1, history2));
		}
	}
}

std::string header()
{
	std::ostringstream s;
	for(auto w : { "time", "sc", "sc2" })
		for(auto m : { "hb", "heur", "syn", "synx" })
			s << "\t" << w << "_" << m << "\tsd_" << w << "_" << m;
	return s.str();
}

void mutation_loop(std::string datadir, double indelrate)
{
	size_t bslen = default_bslen;
	int iters = default_iters;

	std::ostringstream fname;
	fname << datadir << "/mutationrate_i" << indelrate << ".out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#mutrate" << header() << std::endl;

	for(double mut = 0; mut < .31; mut += .02)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, bslen, mut, indelrate, data);
		print(out, mut, data);
	}
}

void indel_loop(std::string datadir, double mutrate)
{
	size_t bslen = default_bslen;
	int iters = default_iters;

	std::ostringstream fname;
	fname << datadir << "/indelrate_m" << mutrate << ".out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#indelrate" << header() << std::endl;

	for(double indel = 0; indel < .01; indel += .0005)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, bslen, mutrate, indel, data);
		print(out, indel, data);
	}
}

void length_loop(std::string datadir, double mutrate, double indelrate)
{
	int iters = 10;

	std::ostringstream fname;
	fname << datadir << "/lengths_m" << mutrate << "_i" << indelrate << ".out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#bslen" << header() << std::endl;

	for(double len = 1e3; len < 1e6; len *= 2.1544)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, len, mutrate, 0, data);
		print(out, len, data);
	}
}

void synapsing_loop(std::string datadir, double mutrate, double indelrate)
{
	int iters = default_iters;

	std::ostringstream fname;
	fname << datadir << "/synapsing.out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#synlen" << header() << std::endl;

	synapsing_alignlen = 20;
	for(synapsing_minlen = 1; synapsing_minlen <= 10; synapsing_minlen++)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, default_bslen, mutrate, indelrate, data, true);
		print(out, synapsing_minlen, data);
	}
}

void synapsing_loop_2(std::string datadir, double mutrate, double indelrate)
{
	int iters = default_iters;

	std::ostringstream fname;
	fname << datadir << "/synapsing_alignlen.out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#synlen" << header() << std::endl;

	synapsing_minlen = 3;
	for(synapsing_alignlen = 1; synapsing_alignlen <= 32; synapsing_alignlen *= 2)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, default_bslen, mutrate, indelrate, data, true);
		print(out, synapsing_alignlen, data);
	}
}

std::vector<std::string> explode(const std::string &str)
{
	std::istringstream ss(str);
	std::vector<std::string> v;
	while(ss)
	{
		std::string s;
		ss >> s;
		if(ss)
			v.push_back(s);
	}
	return v;
}

std::vector<std::string> get_words(const bitstring &bs)
{
	std::vector<std::string> ret;
	for(size_t i = 0; i + 16 <= bs.bits(); ++i)
	{
		char c = (char)bs.get(i, 8);
		if(c != '{')
			continue;
		std::vector<char> s;
		for(i += 8; i + 8 <= bs.bits(); i += 8)
		{
			c = (char)bs.get(i, 8);
			if(c == '}')
				break;
			s.push_back(c);
		}
		ret.push_back(std::string(s.begin(), s.end()));
	}
	return ret;
}

static int bigneg = -1000000;

static inline int letter_distance(char target, char letter)
{
	char c = target ^ letter;
	if(fitness_function == fit_orig)
		return __builtin_popcount(c);

	if(fitness_function == fit_barrier4)
	{
		int dist = 0;
		for(int i = 0; i < 2; ++i)
		{
			char d = c & 0xf;
			if(__builtin_popcount(d) == 1)
				dist += 3;
			else if(!d)
				dist += 1;
			c >>= 4;
		}
		return dist;
	}

	int pc = __builtin_popcount(c);
	if(pc <= 2)
		return pc;
	else if(pc <= 3)
		return 6;
	return 3;
}

int word_distance(const std::string &target, const std::string &text)
{
	int gapcost = 10;
	size_t n = target.size() + 1, m = text.size() + 1;
	static std::vector<int> score;
	score.resize(n * m);
	for(size_t i = 0, ix = 0; i < n; ++i)
	{
		for(size_t j = 0; j < m; ++j, ++ix)
		{
			int ins = i ? score[ix - m] - gapcost : bigneg;
			int del = j ? score[ix - 1] - gapcost : bigneg;
			int match = i && j ? score[ix - m - 1] -
				letter_distance(target[i-1], text[j-1]) :
				(i || j ? bigneg : 0);
			score[ix] = std::max(ins, std::max(del, match));
		}
	}
	return score.back();
}

int text_distance(const std::vector<std::string> &target,
	const std::vector<std::string> &text)
{
	int gapcost = 0;
	int lettercost = 10;
	size_t n = target.size() + 1, m = text.size() + 1;
	static std::vector<int> score;
	score.resize(n * m);
	for(size_t i = 0, ix = 0; i < n; ++i)
	{
		for(size_t j = 0; j < m; ++j, ++ix)
		{
			int ins = i ? score[ix - m] - gapcost - lettercost * target[i-1].size() : bigneg;
			int del = j ? score[ix - 1] - gapcost - lettercost * text[j-1].size(): bigneg;
			int match = i && j ? score[ix - m - 1] +
				word_distance(target[i-1], text[j-1]) :
				(i || j ? bigneg : 0);
			score[ix] = std::max(ins, std::max(del, match));
		}
	}
	return score.back();
}

static int choose_tournament(const std::vector<int> &fitness, bool maximize)
{
	int index1 = gsl_rng_uniform_int(rng, fitness.size());
	int index2 = gsl_rng_uniform_int(rng, fitness.size() - 1);
	if(index2 >= index1)
		++index2;
	if(fitness[index1] == fitness[index2])
		return gsl_rng_uniform_int(rng, 2) ? index1 : index2;
	if((fitness[index1] > fitness[index2]) == maximize)
		return index1;
	return index2;
}

void crossover(const bitstring &src1, const bitstring &src2,
	bitstring &dest1, bitstring &dest2, double crossover_cross_rate,
	int method)
{
	std::vector<char> a1, a2;
	if(!method)
		align_hirschberg(src1, src2, a1, a2);
	else if(method == 1)
		align_heuristic(src1, src2, a1, a2, heuristic_maxdiff);
	else if(method == 2)
		align_synapsing(src1, src2, a1, a2, synapsing_minlen, 0);
	else
		align_synapsing(src1, src2, a1, a2, synapsing_minlen, synapsing_alignlen);

	std::vector<char> b1, b2;

	// Special case: if 50% crossover chance, possibly cross at every
	// mismatch and gap start
	if(crossover_cross_rate == .5)
	{
		bool prevgap = false;
		for(size_t i = 0; i < a1.size(); i++)
		{
			bool curgap = a1[i] == '-' || a2[i] == '-';
			if(((curgap && !prevgap) || (!curgap && a1[i] != a2[i])) &&
				gsl_rng_uniform_int(rng, 2))
			{
				a1.swap(a2);
			}
			prevgap = curgap;
			if(a1[i] != '-')
				b1.push_back(a1[i]);
			if(a2[i] != '-')
				b2.push_back(a2[i]);
		}
	}
	else
	{
		size_t len = gsl_ran_geometric(rng, crossover_cross_rate);
		bool prevalign = false;
		for(size_t i = 0; i < a1.size(); i++)
		{
			bool align = a1[i] != '-' && a2[i] != '-';
			if(align || prevalign)
			{
				if(!--len)
				{
					len = gsl_ran_geometric(rng, crossover_cross_rate);
					a1.swap(a2);
				}
			}
			prevalign = align;
			if(a1[i] != '-')
				b1.push_back(a1[i]);
			if(a2[i] != '-')
				b2.push_back(a2[i]);
		}
	}

	// Convert back to bitstring and return offspring.
	dest1 = bitstring(b1);
	dest2 = bitstring(b2);
}

struct evolver
{
	static double mutation_rate;
	static double indel_prob;
	static double indel_power;
	static double crossover_cross_rate;
	static double crossover_prob;
	static size_t max_length;

	int method;
	std::vector<bitstring> pop;
	std::vector<int> fitness;
	std::vector<std::string> target;
	std::ofstream out;
	double cpu;
	bool done;

	static void setup(double mutrate, double indelp, double xrate, double xprob,
		size_t maxlen)
	{
		mutation_rate = mutrate;
		indel_prob = indelp;
		indel_power = -1.5;
		crossover_cross_rate = xrate;
		crossover_prob = xprob;
		max_length = maxlen;
	}

	evolver(){}

	evolver(const evolver &o):
		method(o.method), pop(o.pop), fitness(o.fitness), target(o.target),
		cpu(0), done(false)
	{
	}

	void init(int method, int bits, int popcnt, const std::string &targetstr)
	{
		this->method = method;
		cpu = 0;
		done = false;

		target = explode(targetstr);
		pop.resize(popcnt, bitstring(bits));
		fitness.resize(popcnt);
		for(int i = 0; i < popcnt; ++i)
		{
			pop[i].rand();
			fitness[i] = text_distance(target, get_words(pop[i]));
		}
	}

	static void log_settings_header(std::ostream &out)
	{
		out << "#mutation_rate\tindel_prob\tindel_power\tcrossover_cross_rate\t"
			"crossover_prob";
	}

	static void log_settings(std::ostream &out)
	{
		out << mutation_rate << "\t" << indel_prob << "\t" <<
			indel_power << "\t" << crossover_cross_rate << "\t" << crossover_prob;
	}

	void open_log(const std::string &fname, bool log_crossovers = false)
	{
		out.open(fname);
		if(!out)
			throw "unable to open " + fname;
		if(!target.empty())
		{
			out << "#target: ";
			std::copy(target.begin(), target.end(),
				std::ostream_iterator<std::string>(out, " "));
			out << "\n";
		}
		log_settings_header(out);
		out << "\n#";
		log_settings(out);
		out << "\n";
		if(log_crossovers)
		{
			out << "#gen\tf_p1\tf_p2\tl_p1\tl_p2"
				"\thi_f_ch1\thi_f_ch2\thi_l_ch1\thi_l_ch2"
				"\theur_f_ch1\theur_f_ch2\theur_l_ch1\theur_l_ch2"
				"\tsyn_f_ch1\tsyn_f_ch2\tsyn_l_ch1\tsyn_l_ch2" <<
				std::endl;
		}
		else
		{
			out << "#generation\tavgfitness\tavgbits\tbestfitness"
				"\tbestbits\tcputime\ttext" << std::endl;
		}
	}

	void step()
	{
		if(gsl_rng_uniform(rng) >= crossover_prob)
		{
			int index = choose_tournament(fitness, true);
			int index_replace;
			while((index_replace = choose_tournament(fitness, false)) == index) ;

			pop[index_replace] = pop[index];
			double m = gsl_rng_uniform(rng);
			if(m < indel_prob * .5 && pop[index_replace].bits() < max_length)
				pop[index_replace].duplication(indel_power, max_length);
			else if(m < indel_prob)
				pop[index_replace].deletion(indel_power, 1);
			else
				pop[index_replace].mutate(mutation_rate);
			
			fitness[index_replace] = text_distance(target, get_words(pop[index_replace]));
		}
		else
		{
			// Draw four different individuals, two good and two bad.
			std::vector<int> ixs;
			while(ixs.size() < 4)
			{
				int ix;
				while(std::find(ixs.begin(), ixs.end(),
					ix = choose_tournament(fitness, ixs.size() < 2)) != ixs.end()) ;
				ixs.push_back(ix);
			}

			// For synapsing, randomize crossover at every point.
//			double xrate = method == 2 ? crossover_cross_rate : .5;
			double xrate = crossover_cross_rate;
			crossover(pop[ixs[0]], pop[ixs[1]], pop[ixs[2]], pop[ixs[3]],
				xrate, method);
			fitness[ixs[2]] = text_distance(target, get_words(pop[ixs[2]]));
			fitness[ixs[3]] = text_distance(target, get_words(pop[ixs[3]]));
		}
	}

	void store_comparison(std::vector<int> &data, int f1, int f2, int l1, int l2)
	{
		if(f1 > f2)
		{
			std::swap(f1, f2);
			std::swap(l1, l2);
		}
		data.push_back(f1);
		data.push_back(f2);
		data.push_back(l1);
		data.push_back(l2);
	}

	void log_bad_crossover(const bitstring &p1, const bitstring &p2,
		const bitstring &c1, const bitstring &c2)
	{
		static int count = 0;
		if(count > 10)
			return;
		++count;
		std::cout << "Bad crossover "
			"(f=" << text_distance(target, get_words(p1)) << ", b=" << p1.bits() << ") + "
			"(f=" << text_distance(target, get_words(p2)) << ", b=" << p2.bits() << ") -> "
			"(f=" << text_distance(target, get_words(c1)) << ", b=" << c1.bits() << ") + "
			"(f=" << text_distance(target, get_words(c2)) << ", b=" << c2.bits() << ")\n";
		p1.print(std::cout);
		std::cout << "\n";
		p2.print(std::cout);
		std::cout << std::endl;
	}

	void compare_crossovers(int gen)
	{
		int ix1 = choose_tournament(fitness, true), ix2;
		while((ix2 = choose_tournament(fitness, true)) == ix1) ;

		std::vector<int> data;
		store_comparison(data, fitness[ix1], fitness[ix2],
			pop[ix1].bits(), pop[ix2].bits());
		bitstring ch1, ch2;
		for(int m = 0; m <= 2; ++m)
		{
			crossover(pop[ix1], pop[ix2], ch1, ch2, crossover_cross_rate, m);
			int cf1 = text_distance(target, get_words(ch1));
			int cf2 = text_distance(target, get_words(ch2));
			store_comparison(data, cf1, cf2, ch1.bits(), ch2.bits());

			if(!method && !m && gen > 20000 && std::max(cf1, cf2) < -200)
				log_bad_crossover(pop[ix1], pop[ix2], ch1, ch2);
		}
		out << gen;
		for(int d : data)
			out << "\t" << d;
		out << std::endl;
	}
	
	void addtime(double t)
	{
		cpu += t;
	}

	bool is_done()
	{
		if(done)
			return true;
		if(*std::max_element(fitness.begin(), fitness.end()) == 0)
			done = true;
		return done;
	}

	const bitstring &get_best() const
	{
		return pop[std::distance(fitness.begin(),
			std::max_element(fitness.begin(), fitness.end()))];
	}

	void log(int generation)
	{
		int best = std::distance(fitness.begin(),
			std::max_element(fitness.begin(), fitness.end()));
		auto words = get_words(pop[best]);
		for(auto &w : words)
		{
			for(auto &c : w)
			{
				if(c < 32 || c >= 127)
					c = '*';
			}
		}

		out << generation << "\t" <<
			 std::accumulate(fitness.begin(), fitness.end(), 0.) /
			 	fitness.size() << "\t" <<
			 std::accumulate(pop.begin(), pop.end(), 0., [](double s,
			 	const bitstring &bs) { return s + bs.bits(); } ) / fitness.size() <<
			 "\t" << fitness[best] << "\t" << pop[best].bits() <<
			 "\t" << cpu << "\t";
		std::copy(words.begin(), words.end(),
			std::ostream_iterator<std::string>(out, " "));
		out << std::endl;
	}
};

double evolver::mutation_rate;
double evolver::indel_prob;
double evolver::indel_power;
double evolver::crossover_cross_rate;
double evolver::crossover_prob;
size_t evolver::max_length;

void evolve(const std::string &datadir, const std::string &targetstr,
	int bits, int popcnt, int generations)
{
	evolver::setup(0.004, 0.5, 0.5, 0.02, bits * genome_maxlen_factor);

	std::vector<evolver> evolvers(2);
	evolvers[0].init(1, bits, popcnt, targetstr);
	evolvers[0].open_log(datadir + "/evolution_heuristic.out");
	evolvers[1].init(2, bits, popcnt, targetstr);
	evolvers[1].open_log(datadir + "/evolution_synapsing.out");

	for(auto &ev : evolvers)
		ev.log(0);
	int log_every = 1000;
	for(int it = 0; it < generations; it += log_every)
	{
		for(auto &ev : evolvers)
		{
			getcputime();
			for(int i = 0; i < log_every; ++i)
				ev.step();
			ev.addtime(getcputime());
			ev.log(it + log_every);
		}
	}	
}

void evolve_comparison(const std::string &datadir, const std::string &targetstr,
	int bits, int popcnt, int generations)
{
	evolver::setup(0.004, 0.5, 0.5, 0.2, bits * genome_maxlen_factor);

	std::vector<evolver> evolvers(1);
	evolvers[0].init(0, bits, popcnt, targetstr);
	evolvers.push_back(evolvers[0]);
	evolvers.push_back(evolvers[0]);
	evolvers[1].method = 1;
	evolvers[2].method = 2;
	evolvers[0].open_log(datadir + "/cross_compare_m0.out", true);
	evolvers[1].open_log(datadir + "/cross_compare_m1.out", true);
	evolvers[2].open_log(datadir + "/cross_compare_m2.out", true);

	for(auto &ev : evolvers)
		ev.compare_crossovers(0);
	int log_every = 10;
	for(int it = 0; it < generations; it += log_every)
	{
		for(auto &ev : evolvers)
		{
			for(int i = 0; i < log_every; ++i)
				ev.step();
			ev.compare_crossovers(it + log_every);
		}
	}	
}

void print_more_stats(std::ostream &out, std::vector<double> &vec)
{
	std::sort(vec.begin(), vec.end());

	double sum = 0, sum2 = 0, n = vec.size();
	for(auto v : vec)
	{
		sum += v;
		sum2 += v * v;
	}
	double mean = sum / n;
	double sd = std::sqrt(std::max(0., sum2 / n - sqr(mean)) / (n - 1));
	double min = vec.front(), max = vec.back();
	double median = .5 * (vec[(vec.size() - 1) / 2] + vec[vec.size() / 2]);

	out << "\t" << mean << "\t" << sd << "\t" << min << "\t" << median << "\t" << max;
}

void evolve_many(const std::string &datadir, const std::string &targetstr,
	int method, int bits, int popcnt, int repeats, double mutrate,
	double indelprob, double xbitrate, double xprob)
{
	evolver::setup(mutrate, indelprob, xbitrate, xprob,
		bits * genome_maxlen_factor);

	std::vector<evolver> evolvers(1);
	std::vector<std::vector<double>> generations(evolvers.size());
	std::vector<std::vector<double>> cpu(evolvers.size());

	std::string fname = datadir + "/evolution_runs.out";
	std::ofstream out(fname);
	if(!out)
		throw "unable to open " + fname;
	out << "#method\tgenerations\tcpu\tbits" << std::endl;

	fname = datadir + "/evolution_summary.out";
	std::ofstream endout(fname);
	if(!endout)
		throw "unable to open " + fname;
	evolver::log_settings_header(endout);
	endout << "\tmethod";
	for(auto &a : {"gen", "cpu"})
		for(auto &b : {"mean", "sd", "min", "med", "max"})
			endout << "\t" << a << "_" << b;
	endout << std::endl;

	int steps_per_turn = 1000;
	for(int rep = 0; rep < repeats; ++rep)
	{
		evolvers[0].init(method, bits, popcnt, targetstr);

		for(int gen = 0; ; gen += steps_per_turn)
		{
			bool any_left = false;
			for(size_t e = 0; e < evolvers.size(); ++e)
			{
				auto &ev = evolvers[e];
				if(ev.is_done())
					continue;
				getcputime();
				int step = 0;
				for(; step < steps_per_turn && !ev.is_done(); ++step)
					ev.step();
				ev.addtime(getcputime());

				if(ev.is_done())
				{
					generations[e].push_back(gen + step);
					cpu[e].push_back(ev.cpu);
					out << ev.method << "\t" << gen + step <<
						"\t" << ev.cpu << "\t" << ev.get_best().bits() << std::endl;
				}
				else
					any_left = true;
			}
			if(!any_left)
				break;
		}
	}

	for(size_t e = 0; e < evolvers.size(); ++e)
	{
		evolver::log_settings(endout);
		endout << "\t" << evolvers[e].method;
		print_more_stats(endout, generations[e]);
		print_more_stats(endout, cpu[e]);
		endout << "\n";
	}
}


int printusage(char const* const* argv)
{
	std::string a = std::string(argv[0]);
	std::cerr <<
		"usage: " << a << " <datadir> <mode> <params>\n"
		"modes:\n"
		"   mut|muti|mutii|len|leni|lenm|indel|indelm|syn|syn2\n"
		"   test\n"
		"   evo fitfunc targetstr bits pop generations\n"
		"   crosscomp fitfunc targetstr bits pop generations\n"
		"   evomany fitfunc target method bits pop replicates mutrate indelprob xbitrate xprob\n"
		"fitfunc: simple|barrier4|barrier8\n";
	return 1;
}

int main(int argc, const char **argv)
{
	gsl_rng_set(rng, randrandseed());

	if(argc < 3)
		return printusage(argv);

	// directory to save/load from
	std::string datadir(argv[1]);
	// Extract/remember the directory we're running from
	sourcedir(argv[0]);

	std::string mode(argv[2]);
	int args;      // correct length of argc for
	int optargs = 0;   // allow some optional arguments

	if(mode == "mut" || mode == "muti" || mode == "mutii" ||
		mode == "len" || mode == "lenm" || mode == "leni" ||
		mode == "indel" || mode == "indelm" ||
		mode == "test" || mode == "syn" || mode == "syn2")
	{
		args = 0;
	}
	else if(mode == "evo" || mode == "crosscomp")
		args = 5;
	else if(mode == "evomany")
		args = 10;
	else
		return printusage(argv);

	if(argc < args + 3 || argc > args + 3 + optargs)
		return printusage(argv);

	int c = 3;
	if(args > 3)
	{
		std::string func(argv[3]);
		if(func == "simple")
			fitness_function = fit_orig;
		else if(func == "barrier4")
			fitness_function = fit_barrier4;
		else if(func == "barrier8")
			fitness_function = fit_barrier8;
		else
			return printusage(argv);
		++c;
	}

	try
	{
		if(mode == "mut")
			mutation_loop(datadir, 0);
		else if(mode == "muti")
			mutation_loop(datadir, 0.001);
		else if(mode == "mutii")
			mutation_loop(datadir, 0.01);
		else if(mode == "len")
			length_loop(datadir, 0.05, 0.002);
		else if(mode == "lenm")
			length_loop(datadir, 0.05, 0);
		else if(mode == "leni")
			length_loop(datadir, 0, 0.002);
		else if(mode == "indel")
			indel_loop(datadir, 0);
		else if(mode == "indelm")
			indel_loop(datadir, 0.05);
		else if(mode == "syn")
			synapsing_loop(datadir, 0.05, 0.002);
		else if(mode == "syn2")
			synapsing_loop_2(datadir, 0.05, 0.002);
		else if(mode == "test")
			test_indel_history();
		else if(mode == "evo")
			evolve(datadir, argv[c+0], touint(argv[c+1]), touint(argv[c+2]), touint(argv[c+3]));
		else if(mode == "crosscomp")
		{
			evolve_comparison(datadir, argv[c+0], touint(argv[c+1]),
				touint(argv[c+2]), touint(argv[c+3]));
		}
		else if(mode == "evomany")
		{
			evolve_many(datadir, argv[c+0], touint(argv[c+1]), touint(argv[c+2]),
				touint(argv[c+3]), touint(argv[c+4]), todouble(argv[c+5]),
				todouble(argv[c+6]), todouble(argv[c+7]), todouble(argv[c+8]));
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
