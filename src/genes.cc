#include "genes.h"
#include "common.h"  //for nextbinary
#include <map>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <limits>

constexpr int gene::cooperativity_dist;
constexpr int gene::exclusion_dist;
constexpr double gene::rnap_dna_energy_min;
constexpr double gene::rnap_dna_energy_max;
constexpr double gene::rnap_tf_energy;
constexpr double gene::cooperativity;
// The maximum number of TFs bound to a cisregion, based on the
// space occupied by each TF.
const int gene::maxbindingsites =
	(cisregion_size + exclusion_dist) / (bindmotif_size + exclusion_dist);

int gene::promoter_start = 106;
int gene::promoter_end = 150;
bool gene::enable_repressors = false;
int gene::maxhamming_dna = 6;             // default. is overwritten by value in settings file
int gene::maxhamming_interaction = 6;     // default. is overwritten by value in settings file
int gene::maxhamming_startsite = 6;       // default. is overwritten by value in settings file

bitstring gene::startsite_tf  = bitstring("01000010110110011011111101101010");
bitstring gene::startsite_deg = bitstring("10011100110010110101000100101001");

/*
constexpr int gene::cisregion_size;
constexpr int gene::bindmotif_size;
constexpr int gene::interactionregion_size;
constexpr int gene::rnap_dna_energy_size;
constexpr int gene::productionrate_size;
constexpr int gene::lightsensitivity_size;
constexpr int gene::lightactivation_size;

constexpr double gene::productionrate_min;
constexpr double gene::productionrate_max;
constexpr double gene::lightsensitivity_min;
constexpr double gene::lightsensitivity_max;
*/

gene::gene():
	_rnapactivities(2 * maxbindingsites + 3), _genome(gene_size)
{
}

gene::gene(const bitstring &g, bool isregulator): gene()
{
	_genome = g;
	_isregulator = isregulator;
	refresh();
	_level = 0;
}


void gene::randomize()
{
	_genome.rand();
	_isregulator = gsl_rng_uniform_int(rng, 2);

	refresh();

	_level = 0;
}

double gene::calcrnapactivity(int tfs) const
{
	double aff = exp(-_rnap_dna_energy - rnap_tf_energy * tfs);
	return aff / (aff + 1.);
}

// convert a bitstring to a value between min and max
static inline double bits_to_value(bitstring::T value, int bits, double min, double max)
{
	double delta = max - min;
	return min + delta * value / (bitstring::T(-1) >> (bitstring::Tbits - bits));
}

// convert a bitstring to a value between min and max
static inline double bits_to_logvalue(bitstring::T value, int bits, double min, double max)
{
	double delta = std::log(max) - std::log(min);
	return std::exp(std::log(min) + delta * value / (bitstring::T(-1) >> (bitstring::Tbits - bits)));
}

void gene::refresh()
{
	//update parameters
	_rnap_dna_energy = bits_to_value(_genome.get(rnap_dna_energy_pos, rnap_dna_energy_size),
		rnap_dna_energy_size, rnap_dna_energy_min, rnap_dna_energy_max);

	_productionrate = bits_to_logvalue(_genome.get(productionrate_pos, productionrate_size),
		productionrate_size, productionrate_min, productionrate_max);

	_lightsensitivity = bits_to_logvalue(_genome.get(lightsensitivity_pos, lightsensitivity_size),
		lightsensitivity_size, lightsensitivity_min, lightsensitivity_max);

	size_t ones = _genome.substr(lightactivation_pos, lightactivation_size).popcount();
	if(ones <= 2)
		_lightactivation = -1; // dark activated
	else if(ones >= 6)
		_lightactivation = 1;  // light activated
	else
		_lightactivation = 0;  // not light sensitive

	size_t oneiness = _genome.substr(degrader_pos, degrader_size).popcount();
	if(oneiness < 4)
		_isdegrader = false;   // binds to protein and keeps it from degrading
	else
		_isdegrader = true;    // binds to protein and degrades it

	for(int i = -maxbindingsites - 1; i <= maxbindingsites + 1; ++i)
		_rnapactivities[i + maxbindingsites + 1] = calcrnapactivity(i);
}



size_t gene::addtfactorsites(const tfactor &tf, size_t index)
{
	assert(!tf.isregulator());
	bindingsite bs;
	size_t oldsites = _sites.size();

	int maxpos = cisregion_pos + cisregion_size - bindmotif_size;

	for(int pos = cisregion_pos; pos <= maxpos; ++pos)
	{
try{
		int dist = _genome.distance(tf.getbindingmotif(), pos);
		if(dist <= maxhamming_dna)
		{
			bs.pos = pos;
			bs.affinity = tf.dnaaffinity(dist);
			bs.tfindex = index;
			_sites.push_back(bs);
		}
}catch(std::string s){throw "in atfs "+s;}
	}
	return _sites.size() - oldsites;
}

int gene::tfactorhamming(const tfactor &tf) const
{
	assert(!tf.isregulator());
	int mindist = bindmotif_size;
	int maxpos = cisregion_pos + cisregion_size - bindmotif_size;
	for(int pos = cisregion_pos; pos <= maxpos; ++pos)
	{
try{
		int dist = _genome.distance(tf.getbindingmotif(), pos);
		if(dist < mindist)
			mindist = dist;
}catch(std::string s){throw "in tfh "+s;}
	}
	return mindist;
}

// Check _all_ hamming distances to a TF
void gene::tfactorhammingdists(const tfactor &tf, std::vector<int> &hammings) const
{
	assert(!tf.isregulator());

	int maxpos = cisregion_pos + cisregion_size - bindmotif_size;
	for(int pos = cisregion_pos; pos <= maxpos; ++pos)
	{
try{
		int dist = _genome.distance(tf.getbindingmotif(), pos);
		if(dist <= maxhamming_dna)
			hammings.push_back(dist);
}catch(std::string s){throw "in tfhd "+s;}
	}
}



void gene::get_bindingstat(bindingstats &stats, bool sametf) const
{
	// last pos is a flag bit, exclude
	for(size_t i = 0; i < _sites.size() - 1; ++i)
	{
		// check sites in front of site i (note: they are sorted on pos)
		for(size_t j = i+1; j < _sites.size() - 1; ++j)
		{
			if(sametf != (_sites[i].tfindex == _sites[j].tfindex))
				continue;
			//check if overlapping
			if(_sites[i].pos + bindmotif_size > _sites[j].pos)
				stats.overlapping++;
			// check if cooperative binding
			else if(_sites[i].pos + bindmotif_size + cooperativity_dist > _sites[j].pos)
				stats.cooperative++;
			else
				stats.other++;
		}
	}
}

void gene::get_tfstat(std::map<int, size_t> &hist_binding_same,
	std::map<size_t, std::vector<std::pair<int, int>>> &hist_sign,
//	std::map<std::pair<int, int>, int> &hist_sign_same,
	size_t size) const
{
	// size is number of genes in the organism.
	// track how many pos & neg bindings this that TF has to this cisregion
	std::vector<std::pair<int, int>> tfcount(size);

	// last pos is a flag bit, exclude
	for(size_t i = 0; i < _sites.size() - 1; ++i)
	{
		bool neg = in_promoter(_sites[i].pos) || is_repressor(_sites[i].pos);
		size_t tf = _sites[i].tfindex;
		if(neg)
			tfcount[tf].second++;
		else
			tfcount[tf].first++;
	}

	for(size_t tf = 0; tf < size; ++tf)
	{
		auto& count = tfcount[tf];
		hist_binding_same[count.first + count.second]++;
		hist_sign[tf].push_back(count);
	}
}


void gene::updatesites(const std::vector<tfactor> &tfactors)
{
	_sites.clear();
	for(size_t i = 0; i < tfactors.size(); i++)
	{
		if(!tfactors[i].isregulator())
			addtfactorsites(tfactors[i], i);
	}
	sort(_sites.begin(), _sites.end());

	// Add a terminator for convenience
	_sites.push_back(bindingsite(-1));
}

// Updates the list of interactions by matching every degrader to
// every possible site.
void gene::updateinteractions(const std::vector<gene> &proteins)
{
	_interactions.clear();
	for(size_t i = 0; i < proteins.size(); ++i)
	{
		auto &p = proteins[i];
		if(!p.isregulator())
			continue;

		int maxpos = interaction_pos + interactionregion_size - bindmotif_size;
		int mindist = bindmotif_size;
		for(int pos = interaction_pos; pos <= maxpos; ++pos)
		{
try{
			mindist = std::min(mindist,
				_genome.distance(p.getbindingmotif(), pos));
}catch(std::string s){throw "in ui "+s;}
		}
		if(mindist <= maxhamming_interaction)
		{
			interaction ion;
			ion.affinity = p.proteinaffinity(mindist);
			ion.proteinindex = i;
			ion.isdegrader = p.isdegrader();
			_interactions.push_back(ion);
		}
	}
}

void gene::randomize_sites(size_t tfactor_count, int maxhamming)
{
	assert(maxhamming <= maxhamming_dna);
	_sites.clear();
	tfactor tf;
	for(size_t i = 0; i < tfactor_count; ++i)
	{
		// try to add a tf until we succeed.
		do
		{
			tf.randomize();
		}while(tfactorhamming(tf) > maxhamming);

		if(!addtfactorsites(tf, i))
			throw std::string("tfactor oddly not added");
	}
	sort(_sites.begin(), _sites.end());

	// Add a terminator for convenience
	_sites.push_back(bindingsite(-1));
}

// Try add tfactor_count tfs once. Some may not bind, then so be it.
// This avoids bias.
size_t gene::try_binding(size_t tfactor_count,
	std::map<int, size_t> &bindSameHist, int maxhamming)
{
	assert(maxhamming <= maxhamming_dna);
	_sites.clear();
	tfactor tf;
	for(size_t i = 0; i < tfactor_count; ++i)
	{
		tf.randomize();
		int s = 0;					  // to how many sites does it bind?
		if(tfactorhamming(tf) <= maxhamming && !(s = addtfactorsites(tf, i)))
			throw std::string("tfactor oddly not added");
		bindSameHist[s]++;
	}
	sort(_sites.begin(), _sites.end());

	// Add a terminator for convenience
	_sites.push_back(bindingsite(-1));

	return _sites.size() -1;  // How many were actually added
}

void gene::getbindingtfs(std::set<size_t> &tfindices) const
{
	for(size_t i = 0; i < _sites.size() - 1; ++i)
		tfindices.insert(_sites[i].tfindex);
}


void gene::getbindingtfspos(std::map<size_t, size_t> &positions) const
{
	for(size_t i = 0; i < _sites.size() - 1; ++i)
		positions[_sites[i].pos]++;
}


int gene::interaction_sign(size_t tfindex) const
{
	int neg = 0, tot = 0;
	for(auto &site: _sites)
	{
		if(site.pos < 0)  // last element is -1, is a terminator for optimization
			break;
		if(site.tfindex != tfindex)
			continue;
		if(in_promoter(site.pos))
			++neg;
		else if(is_repressor(site.pos))
			++neg;
		++tot;
	}

	if(!tot)
		return 0;
	if(neg == 0)
		return 1;
	if(neg == tot)
		return -1;
	return 0;
/*
	if(neg * 3 >= tot * 2)	// At least 2/3 negative
		return -1;
	return neg * 3 <= tot;	// No more than 1/3 negative
*/
}


double gene::getrnapactivity(const std::vector<double> &tflevels) const
{
	int pos = _sites[0].pos;
	if(pos < 0)
		return _rnapactivities[maxbindingsites + 1];

	// General idea: TFs can bind many sites specified in _sites. Many
	// combinations of sites are unavailable due to exclusion (overlapping
	// sites) but there are still many more micro-states of bound TFs than
	// we can handle within reasonable time. Each state has a statistical
	// weight defined by the binding energy between TFs and DNA, the TF
	// concentration and the cooperativity between TFs within a certain
	// distance.
	//
	// The micro-states can be reduced to macro-states that describe how
	// many TFs are bound without caring about identities and positions,
	// plus one state for blocked RNAP.
	//
	// To find the weights of the macro-states we use a dynamic programming
	// algorithm which goes from one end of the regulatory region and keeps
	// track of the statistical weight given where the previous TF was bound.

	const int cdist = cooperativity_dist + bindmotif_size;
	const int exdist = exclusion_dist + bindmotif_size;

	// These vectors are static as an optimization to avoid reallocations
	// since this function is called in the innermost loop of the ODE solver.

	// The possible numbers of activators, including 0, plus terminator +
	// the possible numbers of repressors, plus terminator
	const int tfstates = 2 * maxbindingsites + 3;
	// Midpoint of the above
	const int midstate = maxbindingsites + 1;
	// The number of different binding positions
	const int maxpos = cisregion_pos + cisregion_size - bindmotif_size;

	// Vector of actually used binding positions, so we avoid indexing by
	// positions that aren't used.
	static std::vector<int> positions(maxpos);
	// The statistical weights of having a certain number of TFs bound
	// at/before a certain position.
	// Indexed by [tfstates * pos + midstate + bound]
	// By bound we mean (activators - repressors).
	// Terminated by -1 at bound higher/lower than 0 for each value of pos.
	static std::vector<double> statweights(maxpos * tfstates);
	// The weight of transcription being blocked with the latest
	// TF bound at a certain position.
	static std::vector<double> blockedweightsvec(maxpos);
	// The statistical weights of having a certain number of TFs bound
	// The sum of the weights before the cooperativity zone. That is, the
	// statistical weight of non-cooperative binding.
	static std::vector<double> weightsumsvec(tfstates);
	// Pointer to the midpoint of the above. Above the midpoint are
	// sums corresponding to statweights_pos, below statweights_neg.
	double *weightsums = &weightsumsvec[midstate];
	// The sum of weights between exclusion distance and cooperativity
	// distance. That is, the statistical weight of having something
	// cooperatively binding to the current position.
	static std::vector<double> coopsumsvec(tfstates);
	double *coopsums = &coopsumsvec[midstate];

	// Pointer to the position and weight of the current binding position
	int *nextposptr = &positions[0];
	double *nextwptr = &statweights[midstate];
	// Points to the most distant previous position that should not yet be
	// considered due to the exclusion distance.
	int *exclposptr = nextposptr;
	double *exclwptr = nextwptr;
	// Points to the most distant previous position that should be subject to
	// the cooperativity factor.
	int *coopposptr = nextposptr;
	double *coopwptr = nextwptr;

	// Corresponding pointers for blockedweights
	double *nextbwptr = &blockedweightsvec[0];
	double *exclbwptr = nextbwptr;
	double *coopbwptr = nextbwptr;
	// Corresponding sums
	double coopsum_b = 0;
	double weightsum_b = 0;

	// 1 + the highest number of bound (activating-repressing) TFs
	int maxb = 1;
	// The lowest number of (activating - repressing)
	int minb = -1;
	// Init statistical weights
	weightsums[-2] = weightsums[-1] = weightsums[1] = 0;
	weightsums[0] = 1;
	coopsums[-2] = coopsums[-1] = coopsums[0] = coopsums[1] = 0;
	// Total affinity of TFs for the current site
	double caffinity = 0.;

	for(const bindingsite *bsp = &_sites[0]; ; bsp++)
	{
		const bindingsite &bs = *bsp;
		// Another TF binding to the same position?
		if(bs.pos == pos)
		{
			caffinity += bs.affinity * tflevels[bs.tfindex];
			continue;
		}

		*nextposptr++ = pos;

		// Advance the exclusion pointer
		while(*exclposptr <= pos - exdist)
		{
			// Add to the statistical weight of cooperative binding as
			// positions slide out of the exclusion zone.
			for(double *f = exclwptr, *t = coopsums; *f >= 0.; ++f, ++t)
				*t += *f;
			for(double *f = exclwptr - 1, *t = coopsums - 1; *f >= 0.; --f, --t)
				*t += *f;
			coopsum_b += *exclbwptr;
			exclwptr += tfstates;
			++exclbwptr;
			++exclposptr;
		}
		// Advance the cooperativity zone pointer
		while(*coopposptr <= pos - cdist)
		{
			// Subtract from the statistical weight of cooperative binding as
			// positions slide out of the cooperativity zone, but instead add
			// to the total weight sum.
			for(double *f = coopwptr, *t = coopsums, *tp = weightsums;
				*f >= 0.; ++f, ++t, ++tp)
			{
				*t -= *f;
				*tp += *f;
			}
			for(double *f = coopwptr - 1, *t = coopsums - 1,
				*tp = weightsums - 1; *f >= 0.; --f, --t, --tp)
			{
				*t -= *f;
				*tp += *f;
			}
			coopsum_b -= *coopbwptr;
			weightsum_b += *coopbwptr;
			coopwptr += tfstates;
			++coopbwptr;
			++coopposptr;
		}

		// Blocking binding in the promoter region?
		if(in_promoter(pos))
		{
			double c = coopsum_b, w = weightsum_b;
			for(int b = minb - 1; b < maxb + 1; b++)
			{
				c += coopsums[b];
				w += weightsums[b];
			}
			*nextbwptr = caffinity * (c * cooperativity + w);
			nextwptr[-1] = -1.;
			nextwptr[0] = -1.;
		}
		else
		{
			// Blocked remains blocked but the weight goes up.
			*nextbwptr = caffinity * (coopsum_b * cooperativity + weightsum_b);

			// Do we need to extend the max number of bound TFs?
			if(coopsums[minb - 1] || weightsums[minb - 1])
			{
				--minb;
				assert(-minb <= maxbindingsites);
				weightsums[minb - 1] = coopsums[minb - 1] = 0.;
			}
			if(coopsums[maxb] || weightsums[maxb])
			{
				++maxb;
				assert(maxb <= maxbindingsites);
				weightsums[maxb] = coopsums[maxb] = 0.;
			}

			// Repressor?
			if(is_repressor(pos))
			{
				// Add a repressor, possibly the first
				for(int b = minb; b < maxb; b++)
				{
					nextwptr[b - 1] = caffinity *
						(coopsums[b] * cooperativity + weightsums[b]);
				}
				nextwptr[minb - 2] = -1.;
				nextwptr[maxb - 1] = -1.;
			}
			else
			{
				// Add an activator, possibly the first
				for(int b = minb; b < maxb; b++)
				{
					nextwptr[b + 1] = caffinity *
						(coopsums[b] * cooperativity + weightsums[b]);
				}
				nextwptr[minb] = -1.;
				nextwptr[maxb + 1] = -1.;
			}
		}
		nextwptr += tfstates;
		++nextbwptr;

		if(bs.pos < 0)
			break;
		pos = bs.pos;
		caffinity = bs.affinity * tflevels[bs.tfindex];
	}
	while(coopwptr < nextwptr)
	{
		for(double *f = coopwptr, *t = weightsums; *f >= 0.; ++f, ++t)
			*t += *f;
		for(double *f = coopwptr - 1, *t = weightsums - 1; *f >= 0.; --f, --t)
			*t += *f;
		coopwptr += tfstates;
		weightsum_b += *coopbwptr++;
	}

	double prob = 0;
	double wsum = weightsum_b;
	for(int b = minb - 1; b <= maxb; b++)
	{
		wsum += weightsums[b];
		prob += _rnapactivities[maxbindingsites + 1 + b] * weightsums[b];
	}

	return prob / wsum;
}


double gene::getrnapactivity_safe(const std::vector<double> &tflevels) const
{
	size_t n = _sites.size() - 1;
	if(!n)
		return calcrnapactivity(0);

	const int cdist = cooperativity_dist + bindmotif_size;
	const int exdist = exclusion_dist + bindmotif_size;

	std::vector<double> weightsumsvec(maxbindingsites * 2 + 1, 0);
	double *weightsums = &weightsumsvec[maxbindingsites];
	weightsums[0] = 1;
	double blockedweightsum = 0;

	std::vector<int> onoff(n, 0);
	while(nextbinary(onoff))
	{
		double weight = 1;
		bool excluded = false;		// impossible due to exclusion?
		int tfs = 0;
		bool blocked = false;

		int prevpos = -cisregion_size;
		for(size_t i = 0; i < n; ++i)
		{
			if(!onoff[i])
				continue;

			int pos = _sites[i].pos;
			if(prevpos >= 0 && pos - prevpos < exdist)
			{
				excluded = true;
				break;
			}

			if(is_repressor(pos))
				--tfs;
			else
				++tfs;

			weight *= _sites[i].affinity * tflevels[_sites[i].tfindex];
			if(prevpos >= 0 && pos - prevpos < cdist)
				weight *= cooperativity;

			prevpos = pos;
			if(in_promoter(pos))
				blocked = true;
		}

		if(!excluded)
		{
			if(blocked)
				blockedweightsum += weight;
			else
				weightsums[tfs] += weight;
		}
	}

	double prob = 0;
	double wsum = blockedweightsum;
	for(int b = -maxbindingsites; b <= maxbindingsites; b++)
	{
		prob += calcrnapactivity(b) * weightsums[b];
		wsum += weightsums[b];
	}
	return prob / wsum;
}

void gene::save(std::ostream &out) const
{
	out << "level: " << _level << std::endl;
}

void gene::load(std::istream &in)
{
	std::string s;
	if(!(in >> s) || s != "level:")
		throw std::string("Error when loading gene level.\t" + s);
	in >> _level;

	refresh();
}

void gene::printbindingsites(std::ostream &out) const
{
	for(auto &site: _sites)
	{
		if(site.pos < 0)
			break;
		bool block = site.pos > promoter_start - bindmotif_size;
		out << "TF " << site.tfindex << " at " << site.pos << ", affinity " <<
			site.affinity << " (" << (block ? "neg" : "pos" ) << ")" << "\n";
	}
}

void gene::printrnapactivities(std::ostream &out) const
{
	out << "TFs\tactivity\n";
	for(size_t i = 0; i < _rnapactivities.size(); ++i)
		out << i << "\t" << _rnapactivities[i] << "\n";
}

bool operator== (const gene &gene1, const gene &gene2)
{
   return gene1.distance(gene2) == 0;
}
bool operator< (const gene &gene1, const gene &gene2)
{
   return gene1._genome.bits() < gene2._genome.bits(); 

}


