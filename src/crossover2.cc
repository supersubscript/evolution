#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "common.h"
#include "bitstring.h"
#include "alignment.h"


size_t default_bslen = 5000;
int default_iters = 400;
static int heuristic_maxdiff = 20;
static int heuristic_minmarg = 4;
static int heuristic_retries = 5;
static size_t synapsing_minlen = 3;
static size_t synapsing_alignlen = 20;
static int genome_maxlen_factor = 1000;

static const size_t new_cost_generations = 30000;

enum crossover_method { cm_hirschberg = 0, cm_heuristic,
	cm_synapsing, cm_synapsing_align, cm_heuristic_new, cm_cutnsplice,
	cm_saga, cm_max, cm_viv };

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
	for(size_t i = 0; i < bs.bits(); ++i)
		history[i] = i;
/*	bool prev = bs[0];
	int j = 0;
	for(size_t i = 0; i < bs.bits(); ++i)
	{
		if(bs[i] != prev)
		{
			prev = bs[i];
			++j;
		}
		history[i] = j;
	}*/
}

void indel_with_history(bitstring &bs, std::vector<int> &history, double indel_length_power)
{
	// duplicate/insert AND delete
	size_t len = ran_discrete_power(rng, indel_length_power, bs.bits());
	size_t from = gsl_rng_uniform_int(rng, bs.bits() - len + 1);
	size_t to = gsl_rng_uniform_int(rng, bs.bits() + 1);
	bs.insert(to, bs.substr(from, len));
	history.reserve(history.size() + len);
//	history.insert(history.begin() + to,
//		history.begin() + from, history.begin() + from + len);
	history.insert(history.begin() + to, len, -1);
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

void crossover_with_history(std::vector<char> &a1, std::vector<char> &a2,
	const std::vector<int> &p1, const std::vector<int> &p2,
   std::vector<int> &c1, std::vector<int> &c2, double crossover_cross_rate)
{
	c1.clear();
	c2.clear();

	size_t ix1 = 0, ix2 = 0;
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
				c1.swap(c2);
			}
		}
		prevalign = align;
		if(a1[i] != '-')
			c1.push_back(p1[ix1++]);
		if(a2[i] != '-')
			c2.push_back(p2[ix2++]);
	}
	assert(ix1 == p1.size());
	assert(ix2 == p2.size());
}

void find_crossover_points(int method,
	const bitstring &bs1, const bitstring &bs2,
   size_t &pos1, size_t &pos2)
{
	if(method == cm_cutnsplice)
	{
		pos1 = gsl_rng_uniform_int(rng, bs1.bits());
		pos2 = gsl_rng_uniform_int(rng, bs2.bits());
	}
	else if(method == cm_saga)
	{
		find_crossover_point_saga(bs1, bs2, pos1, pos2);
	}
	else if(method == cm_viv)
	{
	}
}

void crossover_onepoint_with_history(int method,
	const bitstring &bs1, const bitstring &bs2,
	const std::vector<int> &p1, const std::vector<int> &p2,
   std::vector<int> &c1, std::vector<int> &c2)
{
	size_t pos1, pos2;
	find_crossover_points(method, bs1, bs2, pos1, pos2);
	c1.clear();
	c2.clear();
	c1.insert(c1.end(), p1.begin(), p1.begin() + pos1);
	c1.insert(c1.end(), p2.begin() + pos2, p2.end());
	c2.insert(c2.end(), p2.begin(), p2.begin() + pos2);
	c2.insert(c2.end(), p1.begin() + pos1, p1.end());
}

double score_crossover_history(const std::vector<int> &p2,
	const std::vector<int> &c1, const std::vector<int> &c2)
{
	std::vector<int> seen(p2.size());
	for(int i : p2)
	{
		if(i >= (int)seen.size())
			seen.resize(i + 1);
		if(i >= 0)
			seen[i] = 1;
	}
	for(int i : c1)
	{
		if(i >= 0 && i < (int)seen.size() && seen[i] == 1)
			seen[i] = 2;
	}
	int score = 0, norm = 0;
	for(int i : c2)
	{
		if(i >= 0 && i < (int)seen.size())
		{
			if(seen[i] == 2)
				++score;
			if(seen[i] > 0)
				++norm;
		}
	}
//std::cout <<c1.size()<<"\t"<<c2.size()<<"\t"<<score<<"\t"<<norm<<std::endl;
	return (double)score / norm;
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
	double indel_length_power = -2;

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

void align(int method, const bitstring &bs1, const bitstring &bs2,
	std::vector<char> &a1, std::vector<char> &a2)
{
	if(method == cm_hirschberg)
		align_hirschberg(bs1, bs2, a1, a2);
	else if(method == cm_heuristic)
		align_heuristic(bs1, bs2, a1, a2, heuristic_maxdiff);
	else if(method == cm_synapsing)
		align_synapsing(bs1, bs2, a1, a2, synapsing_minlen, 0);
	else if(method == cm_synapsing_align)
		align_synapsing(bs1, bs2, a1, a2, synapsing_minlen, synapsing_alignlen);
	else if(method == cm_heuristic_new)
	{
		align_heuristic_new(bs1, bs2, a1, a2, heuristic_maxdiff,
			heuristic_minmarg, heuristic_retries);
	}
	else if(method == cm_cutnsplice)
		align_cutnsplice(bs1, bs2, a1, a2);
	else if(method == cm_saga)
		align_saga(bs1, bs2, a1, a2);
	else
		throw std::string("undefined method");
}

// data is time / score / newscore => (hirschberg / heuristic / synapsing => (iteration => value))
void time_methods(int iters, size_t bslen, double mutrate, double indelrate,
	std::vector<std::vector<std::vector<double>>> &data)
{
	int iters_slowmethods = 4;
	double indel_length_power = -2;
	double crossover_cross_rate = .02;

	data.clear();
	data.resize(3, std::vector<std::vector<double>>(cm_max));

	for(int i = 0; i < iters; ++i)
	{
		bitstring bs1(bslen);
		bs1.rand();
		std::vector<int> history1;
		prepare_history(bs1, history1);
		bitstring bs2(bs1);
		std::vector<int> history2(history1);

		// Probabilistic rounding.
		int indels = (int)(indelrate * bslen);
		if(gsl_rng_uniform(rng) < indelrate * bslen - indels)
			++indels;
		double mrate = mutrate / (indels + 1);
		for(int j = 0; j < indels; ++j)
		{
			bs2.mutate(mrate);
			indel_with_history(bs2, history2, indel_length_power);
		}
		bs2.mutate(mrate);

		for(int j = 0; j < cm_max; ++j)
		{
			if((j == cm_hirschberg || j == cm_saga) && i % (iters / iters_slowmethods))
				continue;

			std::vector<char> a1, a2;
			getcputime();
			align(j, bs1, bs2, a1, a2);

			data[0][j].push_back(getcputime());
			data[1][j].push_back(score_alignment(a1, a2));

			double cross_rate = .5 - .5 * exp(-2 * crossover_cross_rate *
				std::min(bs1.bits(), bs2.bits()) / (bs1.bits() + bs2.bits() - a1.size()));

			double cs = 0;
			const int n_cross = j < cm_cutnsplice ? 100 : 1;
			for(int k = 0; k < n_cross; ++k)
			{
				std::vector<int> c1, c2;
				if(j >= cm_cutnsplice)
					crossover_onepoint_with_history(j, bs1, bs2, history1, history2, c1, c2);
				else
					crossover_with_history(a1, a2, history1, history2, c1, c2, cross_rate);
				cs += score_crossover_history(history2, c1, c2);
			}
			data[2][j].push_back(cs / n_cross);
		}
	}
}

std::string header()
{
	std::ostringstream s;
	for(auto w : { "time", "sc", "sc2" })
		for(auto m : { "hb", "heur", "syn", "synx", "he1", "cns", "saga", "viv" })
			s << "\t" << w << "_" << m << "\tsd_" << w << "_" << m;
	return s.str();
}

void mutation_loop(std::string datadir, double relindelrate)
{
	size_t bslen = default_bslen;
	int iters = default_iters;

	std::ostringstream fname;
	fname << datadir << "/mutationrate_ri" << relindelrate << ".out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#mutrate" << header() << std::endl;

	for(double mut = 0; mut < .31; mut += .02)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, bslen, mut, mut * relindelrate, data);
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
	int iters = 400;

	std::ostringstream fname;
	fname << datadir << "/lengths_m" << mutrate << "_i" << indelrate << ".out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();
	out << "#bslen" << header() << std::endl;

	for(double len = 1e3; len < 1e6; len *= 2.1544)
	{
		std::vector<std::vector<std::vector<double>>> data;
		time_methods(iters, len, mutrate, indelrate, data);
		print(out, len, data);
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

/*std::vector<std::string> get_words(const bitstring &bs)
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
}*/


static int choose_tournament(const std::vector<double> &fitness, bool maximize)
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
	if(method >= cm_cutnsplice)
	{
		size_t pos1, pos2;
		find_crossover_points(method, src1, src2, pos1, pos2);

		dest1 = src1.substr(0, pos1);
		dest1.append(src2.substr(pos2, src2.bits() - pos2));
		dest2 = src2.substr(0, pos2);
		dest2.append(src1.substr(pos1, src1.bits() - pos1));
		return;
	}

	std::vector<char> a1, a2;
	align(method, src1, src2, a1, a2);

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

const int coding_bits = 10;
const bitstring gene_start("110011");
const int gene_size = gene_start.bits() + 3 * coding_bits;

void get_phenotype_vector(const bitstring &bs,
	std::vector<double> &phenotype)
{
	double scale = 1 << coding_bits;
	std::vector<std::pair<double, double>> deriv_changes;
	for(size_t p = 0; p + gene_size <= bs.bits(); ++p)
	{
		if(bs.distance(gene_start, p))
			continue;
		p += gene_start.bits();
		double mid = bs.get(p, coding_bits) / (scale - 1);
		p += coding_bits;
		double height = (2 * bs.get(p, coding_bits) - (scale - 1)) / (scale - 1);
		p += coding_bits;
		double width = .5 * (bs.get(p, coding_bits) + 1) / scale;
		p += coding_bits - 1;
		double deriv = height / width;
		deriv_changes.push_back(std::make_pair(mid - width, deriv));
		deriv_changes.push_back(std::make_pair(mid, -2 * deriv));
		deriv_changes.push_back(std::make_pair(mid + width, deriv));
	}

	deriv_changes.push_back(std::make_pair(0., 0.));
	deriv_changes.push_back(std::make_pair(1., 0.));
	std::sort(deriv_changes.begin(), deriv_changes.end());

	phenotype.clear();
	double y = 0, deriv = 0, x = 0;
	const double sampstep = 1 / scale;
	double sampx = 0;
	for(auto &c : deriv_changes)
	{
		double dx = c.first - x;

		while(dx && sampx <= c.first && sampx <= 1)
		{
			double sampy = y + deriv * (sampx - x);
			phenotype.push_back(sampy);
			sampx += sampstep;
		}
		if(sampx > 1)
			break;

		y += deriv * dx;
		deriv += c.second;
		x = c.first;
	}
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
	std::vector<double> fitness;
	std::ofstream out;
	double cpu;
	size_t generation, n_cost;
	bool done;

	static void setup(double mutrate, double indelp, double xrate, double xprob,
		size_t maxlen)
	{
		mutation_rate = mutrate;
		indel_prob = indelp;
		indel_power = -2;
		crossover_cross_rate = xrate;
		crossover_prob = xprob;
		max_length = maxlen;
	}

	evolver(){}

	evolver(const evolver &o):
		method(o.method), pop(o.pop), fitness(o.fitness),
		cpu(0), generation(0), done(false)
	{
	}

	void init(int method, int bits, int popcnt)
	{
		this->method = method;
		cpu = 0;
		generation = 0;
		n_cost = 0;
		done = false;

		bitstring b(bits);
		b.rand();
		pop.clear();
		pop.resize(popcnt, b);
		fitness.resize(popcnt);
		for(int i = 0; i < popcnt; ++i)
		{
			pop[i].mutate(mutation_rate);
//			pop[i].rand();
			fitness[i] = calc_fitness(pop[i]);
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
				"\tbestbits\tcputime\tgenes" << std::endl;
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
			
			fitness[index_replace] = calc_fitness(pop[index_replace]);
		}
		else
		{
			// Draw three different individuals, two good and one bad.
			std::vector<int> ixs;
			while(ixs.size() < 3)
			{
				int ix;
				while(std::find(ixs.begin(), ixs.end(),
					ix = choose_tournament(fitness, ixs.size() < 2)) != ixs.end()) ;
				ixs.push_back(ix);
			}

			// For synapsing, randomize crossover at every point.
			double xrate = crossover_cross_rate;
			bitstring tmpchild;
			crossover(pop[ixs[0]], pop[ixs[1]], pop[ixs[2]], tmpchild,
				xrate, method);
			// Pick the only acceptable child or randomize between them
			if(pop[ixs[2]].bits() > max_length ||
				(tmpchild.bits() <= max_length && gsl_rng_uniform_int(rng, 2)))
			{
				pop[ixs[2]] = tmpchild;
			}
			fitness[ixs[2]] = calc_fitness(pop[ixs[2]]);
		}
		++generation;
	}

	void next_costfunction()
	{
		++n_cost;
		for(size_t i = 0; i < pop.size(); ++i)
			fitness[i] = calc_fitness(pop[i]);
	}

	void step(int steps)
	{
		getcputime();
		for(int i = 0; i < steps; ++i)
			step();
		addtime(getcputime());
	}


	// x between 0 and 1
	double target_func(double x)
	{
//		size_t n_cost = generation / new_cost_generations;
//		const double per = 0.0377;
//		double cycles = 6 + 2 * std::abs(2 / per * fmod(n_cost + .75 * per, per) - 1);
		double cycles = 6 + .01 * n_cost;
		return std::sin(2. * M_PI * cycles * x);
	}

	double calc_fitness(const bitstring &bs)
	{
		double err = 0;
for(int j=0;j<100;++j){
		std::vector<double> phenotype;
		get_phenotype_vector(bs, phenotype);

//		double err = 0;
err = 0;
		for(size_t i = 0; i < phenotype.size(); ++i)
			err += sqr(target_func((double)i / (phenotype.size() - 1)) - phenotype[i]);
}
		return -err;
	}

	int gene_count(const bitstring &bs)
	{
		int count = 0;
		for(size_t p = 0; p + gene_size <= bs.bits(); )
		{
			if(!bs.distance(gene_start, p))
			{
				p += gene_size;
				++count;
			}
			else
				++p;
		}
		return count;
	}

	void print_phenotype(std::ostream &out, const bitstring &bs)
	{
		std::vector<double> phenotype;
		get_phenotype_vector(bs, phenotype);

		for(size_t i = 0; i < phenotype.size(); ++i)
		{
			double x = (double)i / phenotype.size();
			out << x << "\t" << phenotype[i] << "\t" << target_func(x) << "\n";
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

/*
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
	}*/

	void compare_crossovers()
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
			int cf1 = calc_fitness(ch1);
			int cf2 = calc_fitness(ch2);
			store_comparison(data, cf1, cf2, ch1.bits(), ch2.bits());

//			if(!method && !m && generation > 20000 && std::max(cf1, cf2) < -200)
//				log_bad_crossover(pop[ix1], pop[ix2], ch1, ch2);
		}
		out << generation;
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

	double best_fit() const
	{
		return *std::max_element(fitness.begin(), fitness.end());
	}

	void log()
	{
		int best = std::distance(fitness.begin(),
			std::max_element(fitness.begin(), fitness.end()));

		out << generation << "\t" <<
			 std::accumulate(fitness.begin(), fitness.end(), 0.) /
			 	fitness.size() << "\t" <<
			 std::accumulate(pop.begin(), pop.end(), 0., [](double s,
			 	const bitstring &bs) { return s + bs.bits(); } ) / fitness.size() <<
			 "\t" << fitness[best] << "\t" << pop[best].bits() <<
			 "\t" << cpu << "\t" << gene_count(pop[best]);
		out << std::endl;
		cpu = 0;
	}

	void save_fit(const std::string &fname)
	{
		std::ofstream out(fname);
		if(!out)
			throw "failed to open" + fname;
		out << "#x\ty\ty_target\n";
		print_phenotype(out, get_best());
	}

	struct stats
	{
		// avgcost, avgbits, bestfit, bestbits, cpu, bestgenes
		// -> data points
		std::vector<std::vector<double>> values;
		stats(): values(6) {};

		void add(int ix, double v)
		{
			values[ix].push_back(v);
		}

		static void header(std::ostream &out)
		{
			for(auto &a : {"avgfit", "avgbits", "fit", "bits", "cpu", "genes"})
				for(auto &b : {"mean", "sd", "med"})
					out << "\t" << a << "_" << b;
			out << std::endl;
		}

		void print_vec(std::ostream &out, std::vector<double> &vals)
		{
			std::sort(vals.begin(), vals.end());
			double s = 0, s2 = 0, n = vals.size();
			for(auto v : vals)
			{
				s += v;
				s2 += v * v;
			}
			double median = .5 * (vals[(vals.size() - 1) / 2] + vals[vals.size() / 2]);
			double mean = s / n;
			double sd = std::sqrt(std::max(0., s2 / n - sqr(mean)) / (n - 1));
			out << "\t" << mean << "\t" << sd << "\t" << median;
		}

		void print(std::ostream &out)
		{
			for(auto & vals : values)
				print_vec(out, vals);
			out << "\n";
		}

	};

	void addstats(stats &s)
	{
		int best = std::distance(fitness.begin(),
			std::max_element(fitness.begin(), fitness.end()));
		s.add(0, std::accumulate(fitness.begin(), fitness.end(), 0.) / fitness.size());
		s.add(1, std::accumulate(pop.begin(), pop.end(), 0., [](double s,
			const bitstring &bs) { return s + bs.bits(); } ) / fitness.size());
		s.add(2, fitness[best]);
		s.add(3, pop[best].bits());
		s.add(4, cpu);
		s.add(5, gene_count(pop[best]));
	}
};

double evolver::mutation_rate;
double evolver::indel_prob;
double evolver::indel_power;
double evolver::crossover_cross_rate;
double evolver::crossover_prob;
size_t evolver::max_length;

void evolve(const std::string &datadir, int bits, int popcnt, int generations,
	double mutrate, double indelprob, double xbitrate, double xprob)
{
	evolver::setup(mutrate, indelprob, xbitrate, xprob,
		bits * genome_maxlen_factor);

	std::vector<evolver> evolvers(2);
	evolvers[0].init(1, bits, popcnt);
	evolvers[0].open_log(datadir + "/evolution_heuristic.out");
	evolvers[1].init(2, bits, popcnt);
	evolvers[1].open_log(datadir + "/evolution_synapsing.out");

	for(auto &ev : evolvers)
		ev.log();
	int log_every = 1000;
	assert(!(new_cost_generations % log_every));
	while((int)evolvers[0].generation < generations)
	{
		for(auto &ev : evolvers)
		{
			ev.step(log_every);

			if(!(ev.generation % new_cost_generations))
			{
				std::ostringstream ss;
				ss << datadir << "/final_fit_m" << ev.method << ".g" << ev.generation << ".out";
				ev.save_fit(ss.str());
				ev.next_costfunction();
			}

			ev.log();
		}
	}	

	evolvers[0].save_fit(datadir + "/final_fit_m0.out");
	evolvers[1].save_fit(datadir + "/final_fit_m1.out");
}

void evolve_comparison(const std::string &datadir,
	int bits, int popcnt, int generations,
	double mutrate, double indelprob, double xbitrate, double xprob)
{
	evolver::setup(mutrate, indelprob, xbitrate, xprob,
		bits * genome_maxlen_factor);

	std::vector<evolver> evolvers(1);
	evolvers[0].init(0, bits, popcnt);
	evolvers.push_back(evolvers[0]);
	evolvers.push_back(evolvers[0]);
	evolvers[1].method = 1;
	evolvers[2].method = 2;
	evolvers[0].open_log(datadir + "/cross_compare_m0.out", true);
	evolvers[1].open_log(datadir + "/cross_compare_m1.out", true);
	evolvers[2].open_log(datadir + "/cross_compare_m2.out", true);

	for(auto &ev : evolvers)
		ev.compare_crossovers();
	int log_every = 100;
	for(int it = 0; it < generations; it += log_every)
	{
		for(auto &ev : evolvers)
		{
			ev.step(log_every);
			ev.compare_crossovers();
		}
	}	
}

void evolve_many(const std::string &datadir,
	int bits, int popcnt, int method, double mutrate,
	double indelprob, double xbitrate, double xprob)
{
	evolver::setup(mutrate, indelprob, xbitrate, xprob,
		bits * genome_maxlen_factor);

	// Number of target switches
	int n_switches = 1;
	int skip_switches = 0;
	const int steps_per_turn = 500;
	assert(!(new_cost_generations % steps_per_turn));

	std::vector<evolver> evolvers(1);
	size_t statpoints = new_cost_generations / steps_per_turn;
//	std::vector<std::vector<evolver::stats>> stats(evolvers.size(),
//		std::vector<evolver::stats>(statpoints + 1));
	std::vector<std::vector<evolver::stats>> stats(evolvers.size(),
		std::vector<evolver::stats>(3));
	std::vector<std::vector<evolver::stats>> longstats(evolvers.size(),
		std::vector<evolver::stats>(statpoints * n_switches + 1));

	std::ostringstream fname;
	fname << datadir << "/evolution_runs_m" << method << ".out";
	std::ofstream out(fname.str());
	if(!out)
		throw "unable to open " + fname.str();

	int repeats = 400;
	for(int rep = 0; rep < repeats; ++rep)
	{
		evolvers[0].init(method, bits, popcnt);
//		evolvers[1].init(2, bits, popcnt);

		for(size_t e = 0; e < evolvers.size(); ++e)
		{
			auto &ev = evolvers[e];
			ev.addstats(longstats[e][0]);
		}

		for(int sw = 0; sw < n_switches; ++sw)
		{
			if(sw >= skip_switches)
			{
				for(size_t e = 0; e < evolvers.size(); ++e)
					evolvers[e].addstats(stats[e][0]);
			}
			for(size_t g = 0; g < statpoints; ++g)
			{
				for(size_t e = 0; e < evolvers.size(); ++e)
				{
					auto &ev = evolvers[e];

					ev.step(steps_per_turn);
					if(sw >= skip_switches)
						ev.addstats(stats[e][g >= statpoints / 2 ? 2 : 1]);
//						ev.addstats(stats[e][g + 1]);
					ev.addstats(longstats[e][sw * statpoints + g + 1]);
					ev.cpu = 0;
				}
			}
			for(auto &ev : evolvers)
				ev.next_costfunction();
			out << rep << "\t" << sw << std::endl;
		}

		for(size_t e = 0; e < evolvers.size(); ++e)
		{
			std::ostringstream fname;
			fname << datadir << "/evolution_folded_m" << evolvers[e].method << ".out";
			std::ofstream endout(fname.str());
			if(!endout)
				throw "unable to open " + fname.str();
			evolver::log_settings_header(endout);
			endout << "\n#";
			evolver::log_settings(endout);
			endout << "\n" << "#gen";
			evolver::stats::header(endout);

			for(size_t i = 0; i < stats[e].size(); ++i)
			{
				endout << i * steps_per_turn *statpoints/2;
				stats[e][i].print(endout);
			}
		}

		for(size_t e = 0; e < evolvers.size(); ++e)
		{
			std::ostringstream fname;
			fname << datadir << "/evolution_average_m" << evolvers[e].method << ".out";
			std::ofstream endout(fname.str());
			if(!endout)
				throw "unable to open " + fname.str();
			endout << "#gen";
			evolver::stats::header(endout);

			for(size_t i = 0; i < longstats[e].size(); ++i)
			{
				endout << i * steps_per_turn;
				longstats[e][i].print(endout);
			}
		}

	}
}

void evolve_until(const std::string &datadir,
	double endfit, int method, double mutrate,
	double indelprob, double xbitrate, double xprob)
{
	int popcnt = 100;
	int bits = 5000;

	evolver::setup(mutrate, indelprob, xbitrate, xprob,
		bits * genome_maxlen_factor);

	evolver evo;
	evolver::stats stats;
	std::vector<double> generations;

	int repeats = 400;
	for(int rep = 0; rep < repeats; ++rep)
	{
		evo.init(method, bits, popcnt);

		getcputime();
		while(evo.best_fit() < endfit)
			evo.step();
		evo.addtime(getcputime());
		evo.addstats(stats);
		generations.push_back(evo.generation);

		std::ostringstream fname;
		fname << datadir << "/generations_f" << endfit << "_m" << evo.method <<
			"_x" << xprob << ".out";
		std::ofstream endout(fname.str());
		if(!endout)
			throw "unable to open " + fname.str();
		evolver::log_settings_header(endout);
		endout << "\n#";
		evolver::log_settings(endout);
		endout << "\n" << "#xprob\tn\tgen_mean\tgen_sd\tgen_med";
		evolver::stats::header(endout);

		endout << xprob << "\t" << generations.size();
		stats.print_vec(endout, generations);
		stats.print(endout);
	}
}


int printusage(char const* const* argv)
{
	std::string a = std::string(argv[0]);
	std::cerr <<
		"usage: " << a << " <datadir> <mode> <params>\n"
		"modes:\n"
		"   mut|muti|len|leni|lenm|indel|indelm\n"
		"   test\n"
		"   evo bits pop generations mutrate indelprob xbitrate xprob\n"
		"   crosscomp bits pop generations mutrate indelprob xbitrate xprob\n"
		"   evomany bits pop method mutrate indelprob xbitrate xprob\n"
		"   generations bits pop method mutrate indelprob xbitrate xprob\n";
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

	if(mode == "mut" || mode == "muti" ||
		mode == "len" || mode == "lenm" || mode == "leni" ||
		mode == "indel" || mode == "indelm" ||
		mode == "test")
	{
		args = 1;
	}
	else if(mode == "evo" || mode == "crosscomp" || mode == "evomany")
		args = 8;
	else if(mode == "generations")
		args = 7;
	else
		return printusage(argv);

	if(argc < args + 2 || argc > args + 2 + optargs)
		return printusage(argv);

	int c = 2;

	try
	{
		if(mode == "mut")
			mutation_loop(datadir, 0);
		else if(mode == "muti")
			mutation_loop(datadir, 0.05);
		else if(mode == "len")
			length_loop(datadir, 0.05, 0.0025);	// ~62.5 generations
		else if(mode == "lenm")
			length_loop(datadir, 0.05, 0);
		else if(mode == "leni")
			length_loop(datadir, 0, 0.002);
		else if(mode == "indel")
			indel_loop(datadir, 0);
		else if(mode == "indelm")
			indel_loop(datadir, 0.05);
		else if(mode == "test")
			test_indel_history();
		else if(mode == "evo")
		{
			evolve(datadir, touint(argv[c+1]), touint(argv[c+2]),
				touint(argv[c+3]), todouble(argv[c+4]), todouble(argv[c+5]),
				todouble(argv[c+6]), todouble(argv[c+7]));
		}
		else if(mode == "crosscomp")
		{
			evolve_comparison(datadir, touint(argv[c+1]), touint(argv[c+2]),
				touint(argv[c+3]), todouble(argv[c+4]), todouble(argv[c+5]),
				todouble(argv[c+6]), todouble(argv[c+7]));
		}
		else if(mode == "evomany")
		{
			evolve_many(datadir, touint(argv[c+1]), touint(argv[c+2]),
				touint(argv[c+3]), todouble(argv[c+4]),
				todouble(argv[c+5]), todouble(argv[c+6]), todouble(argv[c+7]));
		}
		else if(mode == "generations")
		{
			evolve_until(datadir, todouble(argv[c+1]),
				touint(argv[c+2]), todouble(argv[c+3]),
				todouble(argv[c+4]), todouble(argv[c+5]), todouble(argv[c+6]));
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
