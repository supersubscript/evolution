#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <queue>
#include <iterator>
#include <cstring>
#include <cassert>
#include <gsl/gsl_multifit_nlin.h>
#include <iostream>
#include <fstream>

#include "organism.h"
#include "alignment.h"

using std::swap;
using std::shared_ptr;
using std::make_shared;

// Instantiate defaults. is overwritten by value in settings file
double organism::_cost_penalty_genes = 0;
double organism::_cost_penalty_links = 0;
size_t organism::_genome_length_max = 600*100;
size_t organism::_genome_length_init = 300*100;

double organism::_gene_indel_prob = 0.2;
double organism::_indel_length_power = -1.;
double organism::_gene_mutation_prob = 0.001;
double organism::_gene_crossover_prob = 0.05;
double organism::_gene_crossover_length_prob = 0.001;

int organism::_alignment_optimization_maxdiff = 3;


// for clocks, split day up in this many time windows
size_t clock_window_count = 4;

// Must be set explicitly for somewhere, or else error.
organism::cost_type organism::_cost_type = organism::cost_undefined;

size_t organism::_next_light = 0;
std::vector<std::shared_ptr<lightlevelgetter>> organism::_light_conditions;
std::vector<std::shared_ptr<lightlevelgetter>> organism::_all_light_conditions;


struct ode_parameters
{
	ode_parameters(organism &o, const lightlevelgetter& e, size_t d):
		org(o), input(e), day(d) {}

	organism& org;
	const lightlevelgetter& input;
	size_t day;
};

void organism::set_cost_type(cost_type v)
{
	_cost_type = v;


	if(v == cost_clock_ld)
	{

		_light_conditions.assign(std::initializer_list<shared_ptr<lightlevelgetter>>{
			make_shared<lightlevelgetter_ld>(24, 10),
			make_shared<lightlevelgetter_ld>(24, 14) });
/*			make_shared<lightlevelgetter_ld>(24, 6),
			make_shared<lightlevelgetter_ld>(24, 12),
			make_shared<lightlevelgetter_ld>(24, 18) });*/
		_all_light_conditions = _light_conditions;
	}
	else if(v == cost_clock_hf)
	{
		_all_light_conditions.assign(std::initializer_list<shared_ptr<lightlevelgetter>>{
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/Radhr00.midpt", "hf00_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/Radhr00.midpt", "hf00_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/Radhr98.midpt", "hf98_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/Radhr98.midpt", "hf98_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake92_93.midpt", "hf92_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake92_93.midpt", "hf92_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake94_93.midpt", "hf94_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake94_93.midpt", "hf94_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake95_96_04.midpt", "hf95_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake95_96_04.midpt", "hf95_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake97_96.midpt", "hf97_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake97_96.midpt", "hf97_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake99_93_96.midpt", "hf99_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake99_93_96.midpt", "hf99_25", 0.25),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake02_93_04.midpt", "hf02_200", 2),
			make_shared<lightlevelgetter_hf>(
			sourcedir() + "/harvardforest/fake02_93_04.midpt", "hf02_25", 0.25) });
		swaplight();
	}
}

void organism::swaplight()
{
	if(_cost_type == cost_type::cost_clock_hf)
	{
		size_t a = _next_light;
		size_t b = _next_light + 2;

		_light_conditions.assign(_all_light_conditions.begin() + a,
										 _all_light_conditions.begin() + b);

		_next_light = (b == _all_light_conditions.size()) ? 0 : b;
	}
}

organism::organism():
	_odesolver(0)
{
	//store gsl setup parameters for ODE solver:
	_sim_control.eps_abs = 1e-8;
	_sim_control.eps_rel = 1e-6;
	_sim_control.maxsteps = 1000; // max number of steps in ODE solver
	_cost = -1;
}

organism::organism(std::istream &in): organism()
{
	load(in);
}

organism::organism(size_t initial_genome_length): organism()
{
	_genome.reinitialize(initial_genome_length);
	_genome.rand();
	refresh();
}

void organism::reinitialize(size_t initial_genome_length)
{
	_genome.reinitialize(initial_genome_length);
	_genome.rand();
	refresh();
}

void organism::refresh()
{
	std::vector<gene> old_genes;
	_genes.swap(old_genes);
	assert(_genes.empty());

	const size_t spos = gene::get_startsite_pos();

	//Find start sites of TFs and degraders
	for(int i = 0; i <= (int)_genome.bits() - (int)gene::get_gene_size(); )
	{
		if(_genome.distance(gene::get_startsite_tf(), i + spos) <= gene::get_maxhamming_startsite())
		{
			_genes.push_back(gene(_genome.substr(i, gene::get_gene_size()), false));
			i += gene::get_gene_size();
		}
		else if(_genome.distance(gene::get_startsite_deg(), i + spos) <= gene::get_maxhamming_startsite())
		{
			_genes.push_back(gene(_genome.substr(i, gene::get_gene_size()), true));
			i += gene::get_gene_size();
		}
		else
			++i;
	}

	// Restore concentration level to same as (closest) old gene
	for(size_t i = 0; i < _genes.size(); ++i)
	{
		double best = gene::get_gene_size();
		for(size_t j = 0; j < old_genes.size(); ++j)
		{
			double d = _genes[i].distance(old_genes[j]);
			if(d < best)
			{
				best = d;
				_genes[i].setlevel(old_genes[j].getlevel());
			}
		}
	}

	// Update binding sites
	for(size_t i = 0; i < size(); ++i)
	{
		_genes[i].updatesites(_genes);
		_genes[i].updateinteractions(_genes);
	}

	//initiate the ODE-solver, will be shared by all copies/clones/mutants
	//of this organism.
	if(!_odesolver || _odesolver->dim() != size())
		_odesolver = std::make_shared<simulator>(size());
}

double organism::getcost_clock_ld(std::vector<size_t> &window_genes,
	bool fixed_output_genes)
{
	size_t window_count = clock_window_count;
	// index [k][d][j][i] is: light_condition, day, time window, gene
	std::vector<std::vector<std::vector<std::vector<double>>>> all_expression_windows;
	if(!getexpressionlevels(window_count, all_expression_windows, _light_conditions, true))
		return -1;

	assert(all_expression_windows.size() == _light_conditions.size());
	assert(all_expression_windows[0].size() == 1);
	assert(all_expression_windows[0][0].size() == window_count);
	assert(all_expression_windows[0][0][0].size() == size());

	if(fixed_output_genes)
		assert(window_genes.size() == window_count);
	else
		window_genes.resize(window_count);

	// sum over k to get each gene (i) expression over all the time
	// windows (j) for different light conditions (k).
	std::vector<std::vector<double>> tot_expression_levels(window_count,
		std::vector<double>(size()));
	for(size_t k = 0; k < _light_conditions.size(); ++k)
		for(size_t j = 0; j < window_count; ++j)
			for(size_t i = 0; i < size(); ++i)
				tot_expression_levels[j][i] += std::min(1., all_expression_windows[k][0][j][i]);

	// find max expression in each window:
	std::vector<double> best_cost(window_count, 0);
	for(size_t j = 0; j < window_count; ++j)
	{
		if(fixed_output_genes)
		{
			best_cost[j] = tot_expression_levels[j][window_genes[j]] /
				_light_conditions.size();
		}
		else
		{
			auto me = std::max_element(tot_expression_levels[j].begin(),
				tot_expression_levels[j].end());
			best_cost[j] = *me / _light_conditions.size();
			window_genes[j] = std::distance(tot_expression_levels[j].begin(), me);
		}
	}

	// sum them up!
	double cost = std::accumulate(best_cost.begin(), best_cost.end(), 0.) /
		window_count;

	assert(cost >= 0 && cost <= 1);
	return cost;
}

// Detection of "summer"
static int clock_hf_target(size_t day, size_t days)
{
	double d = double(day) / days + 10. / 365; // Offset for 1 jan vs solstice
	return d > .2 && d < .8;
}


double organism::getcost_clock_hf(std::vector<size_t> &season_genes,
	bool fixed_output_genes)
{
	size_t windows_per_day = 4;

	// index [k][d][j][i] is: light_condition, day, time window, gene
	std::vector<std::vector<std::vector<std::vector<double>>>> all_expression_windows;
	if(!getexpressionlevels(windows_per_day, all_expression_windows, _light_conditions, true))
		return -1;

	assert(all_expression_windows.size() == _light_conditions.size());
	assert(all_expression_windows[0][0].size() == windows_per_day);
	assert(all_expression_windows[0][0][0].size() == size());

	if(fixed_output_genes)
		assert(season_genes.size() == 1);
	else
		season_genes.resize(1);

	// sum over k to get each gene (i) expression over all the days
	// (d) for different light conditions (k).
	std::vector<double> summer_levels(size());

	double tot_days = 0, tot_summer_days = 0;
	for(size_t k = 0; k < _light_conditions.size(); ++k)
	{
		size_t days = _light_conditions[k]->get_days();
		// Find the target expression level for "summer"
		for(size_t d = 0; d < days; ++d)
			tot_summer_days += clock_hf_target(d, days);
		tot_days += days;
	}
	double summer_target = tot_days / tot_summer_days / windows_per_day;

	for(size_t k = 0; k < _light_conditions.size(); ++k)
	{
		size_t days = _light_conditions[k]->get_days();
		for(size_t d = 0; d < days; ++d)
		{
			if(!clock_hf_target(d, days))
				continue;
			for(size_t j = 0; j < windows_per_day; ++j)
			{
				for(size_t i = 0; i < size(); ++i)
				{
					summer_levels[i] +=
						std::min(summer_target, all_expression_windows[k][d][j][i]);
				}
			}
		}
	}

	// find max expression in each window:
	double cost;
	if(fixed_output_genes)
		cost = summer_levels[season_genes[0]] / tot_days;
	else
	{
		auto me = std::max_element(summer_levels.begin(), summer_levels.end());
		cost = *me / tot_days;
		season_genes[0] = std::distance(summer_levels.begin(), me);
	}

	assert(cost >= 0);
	assert(cost <= 1);
	return cost;
}

double organism::getcost(std::vector<size_t> &output_genes, bool fixed_output_genes)
{
	double cost = 0;
	if(_cost_type == cost_clock_ld)
		cost = getcost_clock_ld(output_genes, fixed_output_genes);
	else if(_cost_type == cost_clock_hf)
		cost = getcost_clock_hf(output_genes, fixed_output_genes);
	else
		throw std::string("invalid cost_type in getcost");

	// Simulation failed?
	if(cost < 0)
		return 1;

	// Apply penalties for too many genes or links
	cost /= 1 + _cost_penalty_genes * size();
	if(_cost_penalty_links)
		cost /= 1 + _cost_penalty_links * getexcessivelinkcount();

	// cost goes to zero:
	cost = 1.0 - cost;

	return cost;
}

double organism::getcost(bool forceEval)
{
	//we already know the cost?
	if(_cost != -1 && !forceEval)
		return _cost;

	std::vector<size_t> dummy;
	_cost = getcost(dummy, false);
	return _cost;
}

void organism::find_used_nodes(const std::vector<size_t> &output_genes,
	std::set<size_t> &used_genes) const
{
	std::vector<std::set<size_t>> inputs(size());
	for(size_t i = 0; i < size(); ++i)
	{
		_genes[i].getbindingtfs(inputs[i]);
		_genes[i].getbindingproteins(inputs[i]);
	}
	std::set<size_t> cur_genes(output_genes.begin(), output_genes.end());
	used_genes = cur_genes;
	while(!cur_genes.empty())
	{
		auto g = cur_genes.begin();
		for(auto h : inputs[*g])
		{
			if(used_genes.insert(h).second)
				cur_genes.insert(h);
		}
		cur_genes.erase(g);
	}
}

struct linkcost
{
	size_t from, to;
	bool tf;
	double delta;

	linkcost(size_t f, size_t t, bool tf, double d): from(f), to(t), tf(tf), delta(d) {}

	bool operator<(const linkcost &o) const
	{
		return delta < o.delta || (delta == o.delta &&
			(from < o.from || (from == o.from &&
			(to < o.to || (to == o.to &&
			(tf < o.tf))))));
	}
};

double organism::prune(double cost_increase_limit)
{
	// Find base cost and output genes
	std::vector<size_t> output_genes;
	double curcost = getcost(output_genes, false);
	if(_cost == -1)
		_cost = curcost;

	// Find genes that are relevant to the output
	std::set<size_t> used_genes;
	find_used_nodes(output_genes, used_genes);

	// Find cost delta for all links (and remove unused ones)
	std::vector<linkcost> linkcosts;
	for(size_t i = 0; i < size(); ++i)
	{
		// Unused gene?
		if(used_genes.find(i) == used_genes.end())
		{
			_genes[i].clear_sites();
			_genes[i].clear_interactions();
		}
		else
		{
			std::set<size_t> indices;
			_genes[i].getbindingtfs(indices);
			for(size_t f : indices)
			{
				auto old = _genes[i].remove_sites(f);
				double delta = getcost(output_genes, true) - _cost;
				linkcosts.push_back(linkcost(f, i, true, delta));
				_genes[i].restore_sites(old);
			}
			indices.clear();
			_genes[i].getbindingproteins(indices);
			for(size_t f : indices)
			{
				auto old = _genes[i].remove_interaction(f);
				double delta = getcost(output_genes, true) - _cost;
				linkcosts.push_back(linkcost(f, i, false, delta));
				_genes[i].restore_interactions(old);
			}
		}
	}

	double lastcost = _cost;

	while(true)
	{
		bool altered_anything = false;
		// Sort links and prune one at a time
		std::sort(linkcosts.begin(), linkcosts.end());
		std::vector<linkcost> updatedlinkcosts;

		for(auto &lc : linkcosts)
		{
			if(lc.tf)
			{
				auto old = _genes[lc.to].remove_sites(lc.from);
				double newcost = getcost(output_genes, true);
				if(newcost - _cost > cost_increase_limit)
				{
					_genes[lc.to].restore_sites(old);
					updatedlinkcosts.push_back(
						linkcost(lc.from, lc.to, lc.tf, newcost - _cost));
				}
				else
				{
					lastcost = newcost;
					altered_anything = true;
				}
			}
			else
			{
				auto old = _genes[lc.to].remove_interaction(lc.from);
				double newcost = getcost(output_genes, true);
				if(newcost - _cost > cost_increase_limit)
				{
					_genes[lc.to].restore_interactions(old);
					updatedlinkcosts.push_back(
						linkcost(lc.from, lc.to, lc.tf, newcost - _cost));
				}
				else
				{
					lastcost = newcost;
					altered_anything = true;
				}
			}
		}
		if(!altered_anything)
			break;
		updatedlinkcosts.swap(linkcosts);
	}
	return lastcost;
}

int organism::getlinkcount(bool regulators, bool degraders, const std::set<size_t> *used_genes) const
{
	int links = 0;
	for(size_t i = 0; i < size(); ++i)
	{
		if(used_genes && used_genes->find(i) == used_genes->end())
			continue;
		std::set<size_t> indexes;
		if(regulators)
			_genes[i].getbindingproteins(indexes, degraders);
		else
			_genes[i].getbindingtfs(indexes);
		links += indexes.size();
	}
	return links;
}

int organism::getexcessivelinkcount() const
{
	int links = 0;
	for(auto &c : _genes)
	{
		std::set<size_t> tfindexes;
		c.getbindingtfs(tfindexes);
		if(tfindexes.size() > 2)
			links += std::pow((int)tfindexes.size() - 2, 3);
	}
	return links;
}


// compute cost from two vectors of not neccessarily equal length
static double degree_distance(const std::vector<size_t> &short_vec,
	const std::vector<size_t> &long_vec)
{
	assert(long_vec.size() > 0);
	assert(long_vec.size() >= short_vec.size());
	double sum = 0;

	for(size_t i = 0; i < long_vec.size(); i++)
	{
		for(size_t x = 0; x < short_vec.size(); x++)
		{
			int j = (i * short_vec.size() + x) / long_vec.size();
			sum += sqr((double)long_vec[i] - (double)short_vec[j]);
		}
	}
	double sequence_cost = sum / (1. + short_vec.size() * long_vec.size());
	double length_cost = sqr(long_vec.size() - short_vec.size());

	return sequence_cost + length_cost;
}

double organism::getcost_null(const std::vector<size_t> &m1_in,
	const std::vector<size_t> &m1_out) const
{
	std::vector<size_t> m2_out, m2_in;
	get_degree_distribution(m2_in, m2_out);

	// always have same number of CIS and TFs
	assert(m1_in.size() == m1_out.size());
	assert(m2_in.size() == m2_out.size());

	if(m1_in.size() <= m2_in.size())
		return degree_distance(m1_in, m2_in) + degree_distance(m1_out, m2_out);
	else
		return degree_distance(m2_in, m1_in) + degree_distance(m2_out, m1_out);
}

static int expfunction(const gsl_vector *x, void *params, gsl_vector *f)
{
	const std::vector<double> &data = *(const std::vector<double> *)params;
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	for(size_t i = 0; i < data.size(); ++i)
		gsl_vector_set(f, i, (a * std::exp(b * i)) - data[i]);
	return GSL_SUCCESS;
}

static int expderiv(const gsl_vector *x, __attribute__((unused)) void *params, gsl_matrix *J)
{
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	for(size_t i = 0; i < J->size1; ++i)
	{
		double tmp = std::exp(b * i);
		gsl_matrix_set(J, i, 0, tmp);
		gsl_matrix_set(J, i, 1, a * i * tmp);
	}
	return GSL_SUCCESS;
}

static int exp_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	expfunction(x, params, f);
	return expderiv(x, params, J);
}

static void detrend_exp(std::vector<double> &data)
{
	// Base initial guess on linear fit
	double sumx = 0;
	for(size_t i = 0; i < data.size(); ++i)
		sumx += data[i];

	double parray[] = { sumx / data.size(), 0 };
	gsl_vector_view params = gsl_vector_view_array(parray, 2);

	gsl_multifit_function_fdf fdf;
	fdf.f = &expfunction;
	fdf.df = &expderiv;
	fdf.fdf = &exp_fdf;
	fdf.n = data.size();
	fdf.p = 2;
	fdf.params = &data;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(
		gsl_multifit_fdfsolver_lmsder, data.size(), 2);
	gsl_multifit_fdfsolver_set(solver, &fdf, &params.vector);

	for(int i = 0; i < 100; ++i)
	{
		if(gsl_multifit_fdfsolver_iterate(solver) ||
			gsl_multifit_test_delta(solver->dx, solver->x, 1e-40, 1e-40))
			break;
	}

	double a = gsl_vector_get(solver->x, 0);
	double b = gsl_vector_get(solver->x, 1);

	for(size_t i = 0; i < data.size(); ++i)
		data[i] -= a * std::exp(b * i) - a * std::exp(b * data.size() / 2.);

	gsl_multifit_fdfsolver_free(solver);
}


static double running_average(const std::vector<double> &data, size_t half_window_width)
{
	double totsum = 0;
	for(size_t i = half_window_width; i < data.size() - 1 - half_window_width; ++i)
		totsum += std::abs(data[i]);
	totsum = std::max(totsum, 1e-5 * (data.size() - 1 - half_window_width * 2));

	double deviation = 0;	// from running average
	for(size_t i = 0; i < data.size() - 1 - 2 * half_window_width; ++i)
	{
		double average = 0, wsum = 0;
		for(size_t j = 0; j < 2 * half_window_width + 1; ++j)
		{
			double t = (double(j) - half_window_width) / half_window_width;
//			double w = .25 * (2. - std::abs(t)) * (2. - std::abs(t));
			double w = 1. - std::abs(t);
//			double w = 1.;
			average += w * data[i + j];
			wsum += w;
		}
		average /= wsum;
		deviation += std::abs(average - data[i + half_window_width]);
	}

	return deviation / totsum;
}

void organism::oscillation(std::vector<double> &best_osc_vec, std::vector<size_t> &best_osc_vec_index)
{
	assert(_cost_type == cost_clock_ld || _cost_type == cost_clock_hf);
	assert(_genes.size() > 0);

	best_osc_vec.clear();
	best_osc_vec_index.clear();

	const int settling_time = 72;   // hours in constant conditions
	const int measuring_time = 72;
	const int half_window = 24;
	const int ll_days = (settling_time + measuring_time + half_window + 23) / 24;

	// index [k][d][j][i] is: light, day, time window, gene
	std::vector<std::vector<std::vector<std::vector<double>>>> expression_windows;

	// Because our steps may become smaller than usual
	auto save_steps = _sim_control.maxsteps;
	_sim_control.maxsteps *= 3;
	std::vector<std::shared_ptr<lightlevelgetter>> light_conditions(
		std::initializer_list<shared_ptr<lightlevelgetter>>{
			make_shared<lightlevelgetter_ldll>(24., 1, 12., ll_days, 1000.),
			make_shared<lightlevelgetter_ldll>(24., 1, 12., ll_days, 0.) });

	bool ok = getexpressionlevels(24, expression_windows, light_conditions, true);

	_sim_control.maxsteps = save_steps;

	if(ok)
	{
		for(auto &day_windows : expression_windows)
		{
			assert(day_windows.size() == 1 + ll_days);

			double best_osc = 0, best_osc_old = 0;
			size_t best_index = -1;
			for(size_t i = 0; i < size(); ++i)
			{
				std::vector<double> expression(measuring_time + half_window * 2);
				for(int j = 0; j < (int)expression.size(); ++j)
				{
					int hour = j + settling_time - half_window;
					expression[j] = day_windows[1 + hour / 24][hour % 24][i];
				}

				detrend_exp(expression);
				best_osc = std::max(running_average(expression, half_window), best_osc);
				if(best_osc > best_osc_old)
				{
					best_osc_old = best_osc;
					best_index = i;
				}
			}
			best_osc_vec.push_back(best_osc);
			best_osc_vec_index.push_back(best_index);
		}
	}
	else
	{
		best_osc_vec.assign(2, -1);
		best_osc_vec_index.assign(2, -1);
	}
}


static size_t maxloops = 1000000;

// Depth first search
static void DFS(const std::vector<std::set<size_t>> &net,
	std::vector<int> &explored, size_t s, size_t &loops)
{
	explored[s] = 2;								// mark node as is/being explored now

	for(auto w : net[s])
	{
		if(explored[w] == 2)
			loops++;									// were just here -> found a loop
		else if(explored[w] < 2)				// if w is not previously explored
			DFS(net, explored, w, loops);
		if(loops > maxloops)
			return;
	}

	explored[s] = 1;
}

size_t organism::countloops(bool transcriptional, bool posttranslational) const
{
	// Each node/vector index is a set of indices of target nodes
	std::vector<std::set<size_t>> graph(_genes.size());

	// Load network into graph data structure
	for(size_t i = 0; i < size(); ++i)
	{
		std::set<size_t> indices;
		if(transcriptional)
			_genes[i].getbindingtfs(indices);
		if(posttranslational)
			_genes[i].getbindingproteins(indices);
		for(auto j : indices)
			graph[j].insert(i);
	}

	// 0 = never visited
	// 1 = visited this pass
	// 2 = is being explored right now
	// 3 = dead, never visit again
	std::vector<int> explored(_genes.size(), 0);

	// Count number of loops: all possible (inc. overlapping) paths
	size_t loops = 0;

	for(size_t i = 0; i < graph.size(); ++i)
	{
		if(explored[i] == 0)
		{
			DFS(graph, explored, i, loops);
			if(loops > maxloops)
				break;

			for(auto &j : explored)
				if(j == 1)
					j = 3;
		}
	}
	return loops;
}


// Return a vector where each element is the number of tf binding to a
// cisregion.
void organism::get_degree_distribution(std::vector<size_t> &m_in,
	std::vector<size_t> &m_out) const
{
	m_in.clear();
	m_out.assign(size(), 0);

	for(size_t i = 0; i < _genes.size(); i++)
	{
		std::set<size_t> tf_indices;
		_genes[i].getbindingtfs(tf_indices);

		m_in.push_back(tf_indices.size());

		for(auto &tf_index : tf_indices)
			++m_out[tf_index];
	}

	// ordered vector easier when comparing different networks.
	std::sort(m_in.begin(), m_in.end());
	std::sort(m_out.begin(), m_out.end());
}


bool organism::getexpressionlevels(size_t window_count,
	std::vector<std::vector<std::vector<std::vector<double>>>> &all_expression_windows,
	const std::vector<std::shared_ptr<lightlevelgetter>> &light_conditions, bool normalize)
{
	assert(_cost_type == cost_clock_ld || _cost_type == cost_clock_hf);
	assert(!light_conditions.empty());
	assert(window_count > 0);

	// index [k][d][j][i] is: light_condition, day, time window, gene
	all_expression_windows.resize(light_conditions.size());

	if(_genes.empty())
		return false;

	std::vector<double> genesum(size());
	size_t totdays = 0;	// total number of days in all conditions

	for(size_t k = 0; k < light_conditions.size(); ++k)
	{
		auto &cond = *light_conditions[k];
		if(!entrain(cond))
			return false;

		if(!solve(all_expression_windows[k], window_count, cond, false, false))
			return false;
		assert(all_expression_windows[k].size() == cond.get_days());
		totdays += cond.get_days();

		for(size_t day = 0; day < cond.get_days(); ++day)
		{
			auto &expression_windows = all_expression_windows[k][day];

			assert(expression_windows.size() == window_count);
			assert(expression_windows[0].size() == size());

			//normalize gene expression across the time windows.
			for(size_t j = 0; j < window_count; ++j)
			{
				for(size_t i = 0; i < size(); ++i)
					genesum[i] += std::max(0., expression_windows[j][i]);
			}
		}
	}

	if(normalize)
	{
		// Normalize to average sum 1 in a day
		for(size_t i = 0; i < size(); ++i)
		{
			double norm = genesum[i] > 0 ? totdays / genesum[i] : 0;
			for(size_t k = 0; k < light_conditions.size(); ++k)
			{
				auto &cond = *light_conditions[k];
				for(size_t day = 0; day < cond.get_days(); ++day)
				{
					for(size_t j = 0; j < window_count; ++j)
					{
						all_expression_windows[k][day][j][i] =
							std::max(0., all_expression_windows[k][day][j][i]) * norm;
					}
				}
			}
		}
	}
	return true;
}

// Finds the min/max protein levels and the maximum transcription rates
// in the current light conditions.
void organism::findmaximumexpression(std::vector<double> &maxexpression,
	std::vector<double> &maxrnapactivity, std::vector<double> &maxdecayrate)
{
	// The time resolution for examining these levels
	size_t window_count = 24 * 60;

	std::vector<std::vector<std::vector<double>>> expression;

	maxexpression.assign(size(), 0);
	maxrnapactivity.assign(size(), 0);
	maxdecayrate.assign(size(), 0);

	auto save_steps = _sim_control.maxsteps;
	_sim_control.maxsteps *= 3;
	for(auto &c : _light_conditions)
	{
		if(!entrain(*c))
			throw std::string("Guru meditation 0x01!");
		if(!solve(expression, window_count, *c, true, false))
			throw std::string("Guru meditation 0x02!");
		for(auto &vec : expression)
		{
			assert(vec.size() == 3);
			for(size_t i = 0; i < size(); ++i)
			{
				maxexpression[i] = std::max(vec[0][i], maxexpression[i]);
				maxrnapactivity[i] = std::max(vec[1][i], maxrnapactivity[i]);
				maxdecayrate[i] = std::max(vec[2][i], maxdecayrate[i]);
			}
		}
	}
	_sim_control.maxsteps = save_steps;
}


// Makes non-binary truth tables for the provided gene, using the provided
// max protein levels to define on/off TF levels that go into the rules.
void organism::computerule(const gene &cisreg,
	const std::vector<double> &maxexpression,
	std::vector<double> &RNAPoutput)
{
	// Maximum output level, for normalization
	double maxactivity = 0;

	// Clear the return vector (easier than resizing to the right size)
	RNAPoutput.clear();

	// get which genes (indices) bind to this cis region
	std::set<size_t> tf_index;
	cisreg.getbindingtfs(tf_index);

	// Levels of all TFs
	size_t tfs = tf_index.empty() ? 0 : *(tf_index.rbegin()) + 1;
	std::vector<double> tflevels(tfs);
	std::vector<int> onoff(tf_index.size(), 0);

	// Examine all combinations of high and low levels of the cisregion's
	// input TFs, and normalize the RNAP activity level.
	// We enumerate the rules backwards, so we reverse RNAPoutput afterwards.
	do
	{
		size_t j = 0;
		for(auto ix : tf_index)
		{
			tflevels[ix] = onoff[j] ? 0. : maxexpression[ix];
			++j;
		}

		double RNAPactivity = cisreg.getrnapactivity(tflevels);
		RNAPoutput.push_back(RNAPactivity);

		if(RNAPactivity > maxactivity)
			maxactivity = RNAPactivity;
	} while(nextbinary(onoff));
	reverse(RNAPoutput.begin(), RNAPoutput.end());

	// We normalize expression levels by the highest level actually observed,
	// but only within certain limits as we still want to consider very low
	// levels as off. 0.1 seems like a conservative lower limit.
//	double normalization = std::max(maxactivity, 0.1);
//	for(auto &i : RNAPoutput)
//		i /= normalization;
}


//mutate random equations n times
void organism::mutate()
{
	// Should we try insertion/deletion?
	double indel = gsl_rng_uniform(rng);

	// 1. duplicate gene region and insert
	if(indel < _gene_indel_prob / 2)
	{
		if(_genome.bits() < _genome_length_max)
		{
			_genome.duplication(_indel_length_power, _genome_length_max);
//         size_t len = ran_discrete_power(rng, _indel_length_power,
//				std::min(_genome.bits(), _genome_length_max - _genome.bits()));
//         size_t from = gsl_rng_uniform_int(rng, _genome.bits() - len + 1);
//         size_t to = gsl_rng_uniform_int(rng, _genome.bits() + 1);
//         _genome.insert(to, _genome.substr(from, len));
		}
	}
	// 2.  or gene deletion
	// else?
	else if(indel < _gene_indel_prob && _genome.bits() > 1)
	{
		_genome.deletion(_indel_length_power, 1);
//		size_t len = ran_discrete_power(rng, _indel_length_power, _genome.bits() - 1);
//		size_t start = gsl_rng_uniform_int(rng, _genome.bits() - len);
//		_genome.remove(start, len);
	}

	// Flip bits with probability _gene_mutation_prob
	_genome.mutate(_gene_mutation_prob);
//	for(size_t b = gsl_ran_geometric(rng, _gene_mutation_prob) - 1;
//		 b < _genome.bits(); b += gsl_ran_geometric(rng, _gene_mutation_prob))
//	{
//		_genome.toggle(b);
//	}

	// stored value is now outdated
	_cost = -1;

	// Update the binding site information since it's probably changed
	// and reinitialize the solver if the number of genes has changed
	refresh();
}

// Creates offspring from two parent organisms by crossover of genomes
void organism::crossover(const organism& source1, const organism& source2,
	organism& dest1, organism& dest2)
{
	if(&source1 == &source2)
	{
		dest1 = source1;
		dest2 = source2;
		return;
	}

	std::vector<char> a, b;
	// Compare bitstrings and return aligned character vectors a and b
	align_heuristic(source1._genome, source2._genome, a, b,
		_alignment_optimization_maxdiff);
	assert(a.size() == b.size());

	size_t len = gsl_ran_geometric(rng, _gene_crossover_length_prob);
	std::vector<char> temp1, temp2;

	// Create intervals and crossover. Only count indices not containing
	// gaps in either string when taking interval length. Temporary vectors
	// are for removing unwanted characters (gaps) after crossover.
	// Crossover can happen before and after all non-gaps.
	bool prevalign = false;	// Was previous position a non-gap?
	for(size_t i = 0; i < a.size(); i++)
	{
		bool align = a[i] != '-' && b[i] != '-';
		if(align || prevalign)
		{
			if(!--len)
			{
				len = gsl_ran_geometric(rng, _gene_crossover_length_prob);
				a.swap(b);
			}
		}
		prevalign = align;
		if(a[i] != '-')
			temp1.push_back(a[i]);
		if(b[i] != '-')
			temp2.push_back(b[i]);
	}
	// Convert back to bitstring and return offspring.
	dest1._genome = bitstring(temp1);
	dest2._genome = bitstring(temp2);

	assert(dest1._genome.bits() + dest2._genome.bits() ==
		source1._genome.bits() + source2._genome.bits());

	// Stored cost values are now outdated
	dest1._cost = -1;
	dest2._cost = -1;

	// Update the binding site information since it's probably changed
	// and reinitialize the solver if the number of genes has changed
	dest1.refresh();
	dest2.refresh();
}


// Solve the ODE system, using GNU Scientific lib.
// If findmaxexpression is true, expression_windows contains a
//	vector of three vectors: (maxexpression, maxrnapactivity, maxdecayrate)
// The outermost level in expression_windows is indexed by day.
bool organism::solve(std::vector<std::vector<std::vector<double>>> &expression_windows,
	size_t window_count, const lightlevelgetter& lightinput, bool findmaxexpression,
	bool updateproteinlevels)
{
	const size_t n = size();

	std::vector<double> vars(n);
	for(size_t i = 0; i < n; ++i)
		vars[i] = _genes[i].getlevel();

	//pointer to our variables, this is where it takes the initial values
	_sim_control.dataptr = &vars[0];

	double astepsize = .1;   //adaptive steplength

	ode_parameters parameters(*this, lightinput, 0);
	//takes: (eq.system, jacobian, dimension, void *param). We use the
	//void pointer to be a struct which contains "this" pointer, so it
	//has an object to operate on.
	const gsl_odeiv_system sys = {odefunc, 0, n, &parameters};

	size_t days = lightinput.get_days();
	expression_windows.resize(days);
	for(size_t d = 0; d < days; ++d)
	{
		parameters.day = d;
		auto &day_window = expression_windows[d];
		if(findmaxexpression)
		{
			day_window.resize(4);
			day_window[0].resize(n);
			day_window[1].resize(n);
			day_window[2].resize(n);
		}
		else
			day_window.resize(window_count);

		double time = 0;         //starting time
		size_t steps = 0;        //the number of steps taken.

		//simulate for a while
		for(size_t i = 0; i < window_count; ++i)
		{
			// Get ready to record expression levels.
			auto &expression = day_window[findmaxexpression ? 3 : i];
			// Zero out and make room for the proteins.
			expression.assign(n, 0);

			// End time of simulations
			_sim_control.maxt = lightinput.get_period() * (i + 1) / window_count;

			if(!_odesolver->simulate(_sim_control, sys, steps, astepsize, time,
				expression))
			{
				return false;
			}

			if(findmaxexpression)
			{
				// Window size, needed for getting average expression from integral
				double window_size = lightinput.get_period() / window_count;

				auto &maxexpression = day_window[0];
				for(size_t j = 0; j < n; ++j)
				{
					// Convert integral of expression to average expression
					double e = expression[j] / window_size;
					if(e > maxexpression[j])
						maxexpression[j] = e;
				}
				auto &maxrnapactivity = day_window[1];
				for(size_t j = 0; j < n; ++j)
				{
					double rnapa = _genes[j].getrnapactivity(vars);
					if(rnapa > maxrnapactivity[j])
						maxrnapactivity[j] = rnapa;
				}
				auto &maxdecayrate = day_window[2];
				for(size_t j = 0; j < n; ++j)
				{
					double dec = _genes[j].getdecayrate(vars);
					if(dec > maxdecayrate[j])
						maxdecayrate[j] = dec;
				}
			}
		}

		if(findmaxexpression)
			day_window.resize(3);
	}

	if(updateproteinlevels)
	{
		for(size_t i = 0; i < n; ++i)
			_genes[i].setlevel(std::max(0., vars[i]));
	}
	return true;
}

// Entrain the ODE system to light conditions
bool organism::entrain(const lightlevelgetter& lightinput)
{
	size_t n = size();
	std::vector<double> vars(n);
	for(size_t i = 0; i < n; ++i)
		vars[i] = _genes[i].getlevel();

	ode_parameters parameters(*this, lightinput, 0);

	//pointer to our variables, this is where it takes the initial values
	_sim_control.dataptr = &vars[0];

	//takes: (eq.system, jacobian, dimension, void *param). We use the
	//void pointer to be a struct which contains "this" pointer, so it
	//has an object to operate on.
	const gsl_odeiv_system sys = {odefunc, 0, n, &parameters};

	// simulator::simulate needs this but we don't use it here.
	std::vector<double> dummy;

	double astepsize = .1;   //adaptive steplength

	//simulate for a while
	for(size_t day = 0; day < 30; ++day)
	{
		// Start and end time of simulations
		double time = 0;
		_sim_control.maxt = parameters.input.get_period();

		//the number of steps taken. Could be useful to know afterwards
		size_t steps = 0;
		if(!_odesolver->simulate(_sim_control, sys, steps, astepsize, time,
			dummy))
		{
			return false;
		}

		double rdiff = 0;
		for(size_t i = 0; i < n; ++i)
		{
			double v1 = _genes[i].getlevel(), v2 = vars[i];
			rdiff += fabs(v1 - v2) / std::max(1e-6, v1 + v2);
			_genes[i].setlevel(v2);
		}

		if(rdiff < 1e-6 * n)
			break;
	}

	for(size_t i = 0; i < n; ++i)
		_genes[i].setlevel(vars[i]);
	return true;
}


//We need a static function in order to call compute(), since it's a member
//function, but since "func" is static it behaves almost like a global
//function, which is what we need, to play well with GSL's C-style.
int organism::odefunc(double t, const double y[], double f[], void *params)
{
	//cast void pointer to struct, to get light, and this pointer to
	//organism i.e. ourselves
	ode_parameters* parameters = (ode_parameters *) params;
	double light = parameters->input.get_light(parameters->day, t);
	return parameters->org.compute(light, y, f);
}


//ode-system, on GSL-friendly form:
int organism::compute(double light, const double y[], double f[])
{
	//build a vector with both input variables, and equation variables
	static std::vector<double> all_variables;

	size_t n = size();
	all_variables.resize(n);
	//get the new variable values from GSL (form y[]) first:
	for (size_t i = 0; i < n; i++)
	{
		//don't allow neg. values
		all_variables[i] = std::max(0.0, y[i]) * _genes[i].getactivity(light);
	}

	for(size_t i = 0; i < n; ++i)
	{
		f[i] = _genes[i].getproductionrate() *
			_genes[i].getrnapactivity(all_variables)
			- y[i] * _genes[i].getdecayrate(all_variables);
	}

	return GSL_SUCCESS;
}


// investigate distribution of hamming distances
void organism::gethammingdistances(std::map<int, size_t> &hamminghistogram) const
{
	for(auto &cisreg : _genes)
	{
		std::set<size_t> tf_index;
		cisreg.getbindingtfs(tf_index);

		for(auto &tf : tf_index)
		{
			std::vector<int> hammings;
			cisreg.tfactorhammingdists(_genes[tf], hammings);

			for(size_t i = 0; i < hammings.size(); i++)
				hamminghistogram[hammings[i]]++;
		}
	}
}

// print statistics for number or nodes, edges, and in/out degree
void organism::get_degreestat(std::ostream &out) const {

	//degree_distribution
	std::vector<size_t> m_in;
	std::vector<size_t> m_out;

	get_degree_distribution(m_in, m_out);

	// all edges going out, must also come back
	size_t sum_in = 0, sum_out = 0;
	for(int n : m_in)  sum_in  += n;
	for(int n : m_out) sum_out += n;

	out << "Nodes: " << _genes.size() << "\n";
	out << "Sum in & out: " << sum_in <<"\t" << sum_out << "\n";
	out << "In degree:\n";
	for(size_t i = 0; i < m_in.size(); ++i)
		out << m_in[i] << " ";
	out << "\n";

	out << "Out degree:\n";
	for(size_t i = 0; i < m_out.size(); ++i)
		out << m_out[i] << " ";
	out << "\n" << std::endl;
}


// investigate frequency of occurrence of different binding interactions
void organism::get_tfstats(std::map<int, size_t> &countSameBinding,
	std::map<std::pair<int, int>, int> &countSameBindingSign,
	std::map<std::pair<int, int>, int> &countGoodBindingSign,
	std::map<size_t, std::vector<size_t>> &poshist,
	bindingstats &samebinding, bindingstats &otherbinding) const
{
	// map from tfindex to pair of pos/neg bindings
	std::map<size_t, std::vector<std::pair<int, int>>> countPosNeg;

	size_t first_gene = 0;

	for(size_t i = first_gene; i < size(); ++i)
	{
		auto &gene = _genes[i];
		gene.get_tfstat(countSameBinding, countPosNeg, size());
		gene.get_bindingstat(samebinding, true);
		gene.get_bindingstat(otherbinding, false);
	}

	// map from sum of pos+neg to different ratios of these
	//std::map<size_t, std::vector<size_t>> poshist;

	for(auto& posneg : countPosNeg)
	{
		size_t pos = 0, neg = 0;
		for(auto &pn : posneg.second)
		{
			countSameBindingSign[pn]++;
			pos += pn.first;
			neg += pn.second;
		}
		size_t sum = pos + neg;
		if(sum > 10 && abs(pos - neg) <= .1 * sum)
		{
			for(auto &pn : posneg.second)
				countGoodBindingSign[pn]++;
		}
		if(poshist[sum].empty())
			poshist[sum].resize(sum+1);  //make room for 0 as well
		poshist[sum][pos]++;
	}
}


void organism::get_bindingposstat(std::map<size_t, size_t> &positions)
{
	// get starting position of each bind site for all cis regions
	for(auto& gene : _genes)
		gene.getbindingtfspos(positions);
}

//save state to out-stream.
void organism::load(std::istream &in)
{
	assert(_cost_type != cost_undefined);

	std::string s;

	if(!(in >> s) || s != "cost:")
		throw std::string("Error when loading fitness.\t" + s);
	in >> _cost;

	if(!(in >> s) || s != "next_light:")
			throw std::string("Error when loading next light.\t" + s);
	in >> _next_light;

	size_t genes_size;
	if(!(in >> s) || s != "number_of_genes:")
		throw std::string("Error when loading number of genes.\t" + s);
	in >> genes_size;

	_genes.resize(genes_size);
	for(size_t i = 0; i < genes_size; ++i)
		_genes[i].load(in);

	if(!(in >> s) || s != "genome:")
		throw std::string("Error when loading genome.\t" + s);
	in >> _genome;

	refresh();
}

void organism::save(std::ostream &out) const
{
	out << "cost: "            << _cost << "\n";
	out << "next_light: "      << _next_light << "\n";

	out << "number_of_genes: " << _genes.size() << "\n";
	for(size_t i = 0; i < _genes.size(); ++i)
		_genes[i].save(out);

	out << "genome:\n"         << _genome << "\n";
	out << std::endl;
}

void organism::savenetwork(std::string savefile, std::string format) const
{
	std::ofstream out(savefile);
	if(format == "csv")
	{
		// CSV format: source_node_id;target_node_id
		for(size_t i = 0; i < _genes.size(); i++)
		{
			std::set<size_t> tf_indices;
			_genes[i].getbindingtfs(tf_indices);

			// cis region i is producing the tf with matching index.
			for(auto &tf_index : tf_indices)
				out << tf_index << ";" << i << std::endl;

			// print nodes with no bindings as well, common in random nets
			if(tf_indices.size() == 0)
				out << i << ";" << std::endl;
		}
	}
	else if(format == "tlp")
	{
		out << "(tlp \"2.0\"\n";
		// TLP format for Tulip: (node <node_id> <node_id> ...)
		out << "(nodes";
		for(size_t i = 0; i < _genes.size(); i++)
			out << " " << i;
		out << ")" << std::endl;

		//TLP format: (edge <edge_id> <source_node_id> <target_node_id>)
		size_t edge_id = 0;
		for(size_t i = 0; i < _genes.size(); i++)
		{
			std::set<size_t> tf_indices;
			_genes[i].getbindingtfs(tf_indices);
			for(auto &tf_index : tf_indices)
			{
				out << "(edge " << edge_id << " " << tf_index
					 << " " << i << ")" << std::endl;
				edge_id++;
			}
		}
		out << ")" << std::endl;
	}
	else if(format == "sif")
	{
		// Simple Interaction Format format For Cytoscape:
		//    <source_node_id> <relasionship type> <target_node_id>
		for(size_t i = 0; i < _genes.size(); i++)
		{
			std::set<size_t> tf_indices;
			_genes[i].getbindingtfs(tf_indices);

			// cis region i is producing the tf with matching index.
			for(auto &tf_index : tf_indices)
				out << tf_index << " " << _genes[i].interaction_sign(tf_index)
					 << " " << i << std::endl;
			// print nodes with no bindings as well, common in random nets
			if(tf_indices.size() == 0)
				out << i << std::endl;
		}
	}
	else if(format == "dot")
	{
		// DOT language for graphviz
		out << "digraph network{\n";
		out << "\tsplines=spline;\n";
		out << "\trankdir=LR;\n";
		out << "\tnode [color=black, shape=circle];\n";
		out << "\tedge [color=blue, style=solid];\n";

		// Set up special style for output node(s)
		size_t out_start, out_stop;
		if(_cost_type == cost_clock_ld)
		{
			// outputs are the first n nodes
			out_start = 0;
			out_stop = clock_window_count;
		}
		else if(_cost_type == cost_clock_hf)
		{
			// output is the first node
			out_start = 0;
			out_stop = 1;
		}
		else
			throw std::string("unrecognized cost function in save net");

		//typeset nodes after degraders, and tfs
		for(size_t i = 0; i < _genes.size(); ++i)
		{
			std::string shape = _genes[i].isdegrader() ? "triangle" : "circle";
			out << "\t" << i << " [shape=" << shape;
			std::string fillcolor = "white";
			if(_genes[i].get_lightactivation() == -1)
				fillcolor = "gray";
			else if(_genes[i].get_lightactivation() == 1)
				fillcolor = "yellow";
			else if(_genes[i].get_lightactivation() != 0)
				throw std::string("light activation must be either: -1,0,1.");
			out << ",fillcolor=" << fillcolor << ",style=\"filled\"";
			if(out_start <= i && i < out_stop)
				out << ",penwidth=5";
			out << "]\n";
		}

		// print transcription network
		for(size_t i = 0; i < _genes.size(); i++)
		{
			std::set<size_t> tf_indices;
			_genes[i].getbindingtfs(tf_indices);

			// cis region i is producing the tf with matching index.
			for(auto &tf_index : tf_indices)
			{
				std::string arrowhead;
				int sign = _genes[i].interaction_sign(tf_index);
				if(sign == -1)     arrowhead = "arrowhead=tee, color=\"red\"";     // negative
				else if(sign == 0) arrowhead = "arrowhead=dot, color=\"gray\"";    // ambigous
				else if(sign == 1) arrowhead = "arrowhead=normal, color=\"blue\""; // positive
				else throw std::string("invalid interaction sign: " + std::to_string(sign));
				out << "\t" << tf_index << " -> " << i << " [" << arrowhead << "];\n";
			}
			// print nodes with no bindings as well, common in random nets
			if(tf_indices.empty())
				out << "\t" << i << ";" << std::endl;
		}

		// print protein interaction network
		for(size_t i = 0; i < _genes.size(); i++)
		{
			std::set<size_t> prot_indices;
			_genes[i].getbindingproteins(prot_indices);

			// cis region i is producing the tf with matching index.
			for(auto &prot_index : prot_indices)
			{
				std::string str = "arrowhead=normal, color=\"orange\", style=dashed";
				out << "\t" << prot_index << " -> " << i << " [" << str;
				out << "];\n";
			}
			// print nodes with no bindings as well,
			if(prot_indices.empty())
				out << "\t" << i << ";" << std::endl;
		}

		out << "}" << std::endl;
	}

	else
		throw "invalid format: " + format;
	out.close();
}

void organism::printtimecourse(const std::string &dirname)
{
	assert(_cost_type == cost_clock_ld || _cost_type == cost_clock_hf);

	// Get index to best oscillating genes
	std::vector<double> tmp;
	std::vector<size_t> best_osc_index;
	oscillation(tmp, best_osc_index);

	// Get output genes
	std::vector<size_t> output_genes;
	getcost(output_genes, false);

	// Add oscillating gene index to vector of best gene indices
	output_genes.insert(output_genes.end(), best_osc_index.begin(), best_osc_index.end());

	std::vector<std::shared_ptr<lightlevelgetter>> light_conditions = _all_light_conditions;

	for(int q = 0; q < 2; ++q)
	{
//		window_count = 4;
	size_t window_count = _cost_type == cost_clock_ld ? 24*12 : 4;

		if(q == 1)
		{
			window_count = 24*60/10;
			light_conditions.clear();
			light_conditions.push_back(
				make_shared<lightlevelgetter_ldll>(24., 1, 12., 7, 1000.));
			light_conditions.push_back(
				make_shared<lightlevelgetter_ldll>(24., 1, 12., 7, 0.));
			light_conditions.push_back(
				make_shared<lightlevelgetter_ldll>(24., 1, 12., 7, 100.));
		}

		// index [i][d][j][g] is: light_condition, day, time window, gene
		std::vector<std::vector<std::vector<std::vector<double>>>> all_expression_windows;
		auto save =  _sim_control.maxsteps;
		_sim_control.maxsteps *= 3;
		bool ok = getexpressionlevels(window_count, all_expression_windows, light_conditions, false);
		_sim_control.maxsteps = save;
		if(!ok)
			throw "failed to integrate system " + dirname;

		std::vector<double> decayrates(size());

		// print to file:
		for(size_t i = 0; i < all_expression_windows.size(); ++i)
		{
			const lightlevelgetter &llg = *light_conditions[i];
			std::ostringstream filename;
			filename << dirname << "/tc." << llg.get_name();
			std::ofstream out(filename.str());
			if(!out)
				throw "failed to open " + filename.str() + " for writing";

			std::ostringstream filename_deg;
			filename_deg << filename.str() << "_deg";
			std::ofstream out_deg(filename_deg.str());
			if(!out_deg)
				throw "failed to open " + filename_deg.str() + " for writing";

			for(size_t d = 0; d < all_expression_windows[i].size(); ++d)
			{
				for(size_t j = 0; j < window_count; ++j)
				{
					double tmid = llg.get_period() * (d + (j + .5) / window_count);
					out << tmid;
					for(auto g : output_genes)
						out << "\t" << all_expression_windows[i][d][j][g];
					out << "\n";

					double light = llg.get_light(d, tmid);
					for(size_t g = 0; g < size(); ++g)
					{
						all_expression_windows[i][d][j][g] *= _genes[g].getactivity(light);
					}

					out_deg << tmid;
					for(size_t g = 0; g < size(); g++)
					{
						out_deg << "\t" << _genes[g].getdecayrate(all_expression_windows[i][d][j]);
					}
					out_deg << "\n";
				}
			}
		}
	}

	std::ostringstream filename;
	filename << dirname << "/lightsensitivity";
	std::ofstream out(filename.str());
	if(!out)
		throw "failed to open " + filename.str() + " for writing";
	out << "# sens\tact\tprodrate\n";
	for(auto &g : _genes)
		out << g.get_lightsensitivity() << "\t" << g.get_lightactivation() << "\t" << g.getproductionrate() << "\n";
}

void organism::printbehaviour(const std::string &dirname)
{
	if(_cost_type == cost_clock_ld || _cost_type == cost_clock_hf)
	{
		printtimecourse(dirname);
	}
}

template<class T>
static void load_thing(std::istream &in, const std::string &name, T &var)
{
	std::string s;
	if(!(in >> s) || s != name || !(in >> var))
		throw "Error when loading " + name + "\t" + s;
}

// CT got lazy...
#define LOAD_THING(x) load_thing(in, #x ":", _ ## x);

void organism::load_globals(const std::string &dir)
{
	std::string fname = dir + "/settings";
	std::ifstream in(fname);
	if(!in)
		throw "failed to read organism settings from " + fname;

	unsigned int tmp_cost_type;
	int pstart, pend;
	bool rep;
	int maxh_cis, maxh_int, maxh_start;

	load_thing(in, "cost-type:", tmp_cost_type);
	load_thing(in, "promoter_start:", pstart);
	load_thing(in, "promoter_end:", pend);
	load_thing(in, "enable_repressors:", rep);
	load_thing(in, "maxhamming_dna:", maxh_cis);
	load_thing(in, "maxhamming_interaction:", maxh_int);
	load_thing(in, "maxhamming_startsite:", maxh_start);
	load_thing(in, "cost_penalty_links:", _cost_penalty_links);
	load_thing(in, "cost_penalty_genes:", _cost_penalty_genes);
	LOAD_THING(gene_indel_prob)
	LOAD_THING(indel_length_power)
	load_thing(in, "gene_mutation_prob:", _gene_mutation_prob);
	load_thing(in, "gene_crossover_prob:", _gene_crossover_prob);
	load_thing(in, "gene_crossover_length_prob:", _gene_crossover_length_prob);
	LOAD_THING(alignment_optimization_maxdiff)
	load_thing(in, "genome_length_max:", _genome_length_max);
	load_thing(in, "genome_length_init:", _genome_length_init);

	gene::set_promoter(pstart, pend);
	gene::set_enable_repressors(rep);
	set_cost_type(static_cast<cost_type>(tmp_cost_type));
	gene::set_maxhamming_dna(maxh_cis);
	gene::set_maxhamming_interaction(maxh_int);
	gene::set_maxhamming_startsite(maxh_start);
}

#define SAVE_THING(x) out << "\n" #x ": " << _ ## x;

void organism::save_globals(const std::string &dir)
{
	std::string fname = dir + "/settings";
	std::ofstream out(fname);
	if(!out)
		throw "failed to open " + fname + " for writing";
	out << "cost-type: " << get_cost_type() << "\n";

	out << "promoter_start: " << gene::get_promoter_start();
	out << "\npromoter_end: " << gene::get_promoter_end();
	out << "\nenable_repressors: " << gene::get_enable_repressors();
	out << "\nmaxhamming_dna: " << gene::get_maxhamming_dna();
	out << "\nmaxhamming_interaction: " << gene::get_maxhamming_interaction();
	out << "\nmaxhamming_startsite: " << gene::get_maxhamming_startsite();
	out << "\ncost_penalty_links: " << _cost_penalty_links;
	out << "\ncost_penalty_genes: " << _cost_penalty_genes;
	SAVE_THING(gene_indel_prob)
	SAVE_THING(indel_length_power)
	out << "\ngene_mutation_prob: " << _gene_mutation_prob;
	out << "\ngene_crossover_prob: " << _gene_crossover_prob;
	out << "\ngene_crossover_length_prob: " << _gene_crossover_length_prob;
	SAVE_THING(alignment_optimization_maxdiff)
	out << "\ngenome_length_max: " << _genome_length_max;
	out << "\ngenome_length_init: " << _genome_length_init << "\n";
}

void organism::log_header(std::ostream &log, bool full)
{
	log << "\tcost\tgenes\tgenome\tTFlinks\tdeglinks\tcomplexlinks\tTF-\tTF0\tTF+\tdeg-\tdeg0\tdeg+\tcomp-\tcomp0\tcomp+";
	if(full)
		log << "\tloopsTF\tloopsProt\tloopsAll\toscLL\toscDD";
	log << std::endl;
}

std::vector<int> organism::get_genetypes(const std::set<size_t> *used_genes) const
{
	std::vector<int> genetypes(9);
	for(size_t i = 0; i < size(); ++i)
	{
		if(used_genes && used_genes->find(i) == used_genes->end())
			continue;
		// gt (gene type) is 0,1,2 if protein is TF, regulator (degrader), regualtor (!degrader)
		int gt = _genes[i].isregulator() ? (_genes[i].isdegrader() ? 1 : 2) : 0;
		// gene type light matches columns in log_header
		int gtl = gt * 3 + _genes[i].get_lightactivation() + 1;
		assert(gtl >= 0 && gtl < 9);
		genetypes[gtl]++;
	}
	return genetypes;
}

void organism::log_line(std::ostream &log, bool full)
{
	if(!full)
	{
		log << "\t" << _cost
			 << "\t" << size()
			 << "\t" << _genome.bits()
			 << "\t" << getlinkcount(false, false)  // 2:nd arg. is a stump
			 << "\t" << getlinkcount(true, true)    // regulator, degrader
			 << "\t" << getlinkcount(true, false);  // regulator, not degrader

		auto genetypes = get_genetypes(nullptr);
		for(auto i : genetypes)
			log << "\t" << i;
	}
	else
	{
		// Find base cost and output genes
		std::vector<size_t> output_genes;
		double cost = getcost(output_genes, false);
		// Find genes that are relevant to the output
		std::set<size_t> used_genes;
		find_used_nodes(output_genes, used_genes);

		log << "\t" << cost
			 << "\t" << used_genes.size()
			 << "\t" << _genome.bits()
			 << "\t" << getlinkcount(false, false, &used_genes)
			 << "\t" << getlinkcount(true, true, &used_genes)
			 << "\t" << getlinkcount(true, false, &used_genes);
		auto genetypes = get_genetypes(&used_genes);
		for(auto i : genetypes)
			log << "\t" << i;

		log << "\t" << countloops(true, false)
			<< "\t" << countloops(false, true)
			<< "\t" << countloops(true, true);
		std::vector<double> expression;
		std::vector<size_t> tmp;
		oscillation(expression, tmp);
		for(auto i : expression)
			log << "\t" << i;
	}
	log << std::endl;
}
        
/*
std::vector<size_t> organism::get_diversity()
{
   std::vector<std::vector<gene>> groups;
   // Maximum difference between genes for claiming familiarity 
   int const max_gene_hamming_diff = 100;
   std::vector<bool> found(_genes.size());
//std::cout << _genes.size() << "\n"; 
   for(size_t i = 0; i < _genes.size(); i++)
   {
      if(found[i])
         continue;
      else
      {
         groups.push_back({_genes[i]});
         found[i] = true;
      }
      for(size_t j = i+1; j < _genes.size(); j++)
      {
         if(found[j])
            continue;
// std::cout <<"distance: "<< _genes[i].distance(_genes[j], max_gene_hamming_diff) << "\n";
         
         if (_genes[i].distance(_genes[j], max_gene_hamming_diff) != -1)
         {   
            found[j] = true;
            for(size_t k = 0; k < groups.size(); k++)
            {
               if(std::find(groups[k].begin(), groups[k].end(), _genes[i]) != groups[k].end())
                  groups[k].push_back(_genes[j]);
            }
         }
      }
   }

   std::vector<size_t> numbers;
   
   for(size_t i = 0; i < groups.size(); ++i)
      numbers.push_back(groups[i].size());
   return numbers;
}
*/


