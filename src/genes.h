#ifndef META_GENES_H
#define META_GENES_H

#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <algorithm>

#include "common.h"
#include "bitstring.h"

class gene
{
public:
	// Model parameters

	// Size of the whole regulatory cis region in bits
	static constexpr int cisregion_size = 256;

	// Size of motif that binds to cis/interaction-region if TF/degrader
	static constexpr int bindmotif_size = 32;


private:
	// Size of startsite, that is between cisregion and interaction region
	static constexpr int startsite_size = 32;

	// Size of the whole protein interaction binding region in bits
	static constexpr int interactionregion_size = 128;

	// parameters are characterized by bit strings:
	static constexpr int rnap_dna_energy_size = 16;
	static constexpr int productionrate_size = 16;
	static constexpr int lightsensitivity_size = 16;
	static constexpr int lightactivation_size = 8;
	static constexpr int degrader_size = 8;

	// Starting position on the genome bit string
	static constexpr int cisregion_pos = 0;
	static constexpr int startsite_pos = cisregion_pos + cisregion_size;
	static constexpr int interaction_pos = startsite_pos + startsite_size;
	static constexpr int bindmotif_pos = interaction_pos + interactionregion_size;
	static constexpr int rnap_dna_energy_pos = bindmotif_pos + bindmotif_size;
	static constexpr int productionrate_pos = rnap_dna_energy_pos + rnap_dna_energy_size;
	static constexpr int lightsensitivity_pos = productionrate_pos + productionrate_size;
	static constexpr int lightactivation_pos = lightsensitivity_pos + lightsensitivity_size;
	static constexpr int degrader_pos = lightactivation_pos + lightactivation_size;

public:
	static constexpr int gene_size = degrader_pos + degrader_size;
private:
	// Limits for the energy of RNAP binding to DNA (without TFs)
	static constexpr double rnap_dna_energy_min = -3.;
	static constexpr double rnap_dna_energy_max = 9.;
	// Limits for the protein production rate (hour ^ -1)
	static constexpr double productionrate_min = 1.;
	static constexpr double productionrate_max = 10.;
	// Limits for the protein light activation ((W/m^2)^-1)
	static constexpr double lightsensitivity_min = 1e-4;
	static constexpr double lightsensitivity_max = 1e-1;

	// The difference one TF makes to the free energy of RNAP binding to DNA
	static constexpr double rnap_tf_energy = -3.;
	// The corresponding probability factor
	static constexpr double cooperativity = exp(-rnap_tf_energy);

	// Start and end of the promoter region where TF binding represses RNAP
	// Promoter 44 bits long, centered => middle third is negative (86, 130)
	// Set promoter_start - promoter_end > bindmotif_size (TODO) to disable
	static int promoter_start;
	static int promoter_end;

	// Start site for TFs and degrader protein
	static bitstring startsite_tf;
	static bitstring startsite_deg;

	//keep track on which gene binds to self
	typedef gene tfactor;

	// Whether to enable TFs that repress RNAP binding based on the local
	// cisregion sequence.
	static bool enable_repressors;

	// The maximum distance of two TFs that bind cooperatively
	static constexpr int cooperativity_dist = 10;
	// The maximum distance at which TFs exclude each other from binding
	static constexpr int exclusion_dist = 0;

	// Maximum mismatch between TF and binding site (for binding to occur)
	static int maxhamming_dna;

	// Maximum mismatch between protein binding motif and protein
	// interaction site
	static int maxhamming_interaction;

	// Maximum mismatch for start site binding sequence
	static int maxhamming_startsite;

	// An upper bound for the number of simultaneously occupied binding
	// sites (the number of TFs bound at any one time)
	static const int maxbindingsites;

	// A site in the region where a specific TF can bind. pos is the position
	// of the first bit of the TF sequence in our sequence.
	struct bindingsite
	{
		bindingsite(int p): pos(p) {}
		bindingsite() {}
		int pos;				// Start position of binding along the bitstring
		size_t tfindex;	// Index of the TF in the vector of TFs
		double affinity;	// Binding affinity
		// Comparison operator for sorting according to position
		bool operator<(const bindingsite &c) const
		{
			return pos < c.pos;
		}
	};


public:

	static inline void set_promoter(int pstart, int pend)
	{
		promoter_start = pstart;
		promoter_end = pend;
	}
	static inline void set_enable_repressors(bool r)
	{
		enable_repressors = r;
	}
	static inline int get_promoter_start()
	{
		return promoter_start;
	}
	static inline int get_promoter_end()
	{
		return promoter_end;
	}
	static inline bool get_enable_repressors()
	{
		return enable_repressors;
	}
	static inline int get_maxhamming_dna()
	{
		return maxhamming_dna;
	}
	static inline void set_maxhamming_dna(int m)
	{
		maxhamming_dna = m;
	}
	static inline int get_maxhamming_interaction()
	{
		return maxhamming_interaction;
	}
	static inline void set_maxhamming_interaction(int m)
	{
		maxhamming_interaction = m;
	}
	static inline int get_maxhamming_startsite()
	{
		return maxhamming_startsite;
	}
	static inline void set_maxhamming_startsite(int m)
	{
		maxhamming_startsite = m;
	}

	static inline bitstring get_startsite_tf()
	{
		return startsite_tf;
	}
	static inline bitstring get_startsite_deg()
	{
		return startsite_deg;
	}

	static inline int get_gene_size()
	{
		return gene_size;
	}
	static inline int get_startsite_pos()
	{
		return startsite_pos;
	}

	int distance(gene g) const
	{
		assert(g._genome.bits() == _genome.bits());
		return _genome.distance(g._genome,0);
	}

	int distance(gene g, int max_hamming) const
	{
		assert(g._genome.bits() == _genome.bits());
		return _genome.distance(g._genome,0, max_hamming);
	}
	struct bindingstats
	{
		double overlapping, cooperative, other;
		bindingstats(): overlapping(0), cooperative(0), other(0) {}
	};

	gene(const bitstring &gene, bool isregulator);

	// Doesn't do much
	gene();

	// Randomizes the genome, and all variables
	void randomize();

	// for statistics about binding overlap
	void get_bindingstat(bindingstats &stats, bool sametf) const;

	// for statistics about binding tfs
	void get_tfstat(std::map<int, size_t> &hist_binding_same,
		std::map<size_t, std::vector<std::pair<int, int>>> &hist_sign,
//		std::map<std::pair<int, int>, int> &hist_sign_same,
		size_t size) const;

	// Updates the gene type, production rate and so forth. Must be called
	// for all genes before and updatesites/updateinteractions.
	void updateparameters();

	// Updates the list of binding sites by matching every TF to every
	// possible site.
	void updatesites(const std::vector<tfactor> &tfactors);

	// Updates the list of interactions by matching every degrader to
	// every possible site.
	void updateinteractions(const std::vector<gene> &proteins);

	// Clear the list of tf-binding sites
	void clear_sites()
	{
		_sites.clear();
		// Add the terminator
		_sites.push_back(bindingsite(-1));
	}

	// Clear the list of protein interactions
	void clear_interactions()
	{
		_interactions.clear();
	}

	// Removes all binding sites bound by the given tf, and returns the old vector of sites
	std::vector<bindingsite> remove_sites(size_t index)
	{
		std::vector<bindingsite> old_sites(_sites);
		_sites.erase(std::remove_if(_sites.begin(), _sites.end() - 1,
			[index](bindingsite &x){return index == x.tfindex;}), _sites.end() - 1);
		return old_sites;
	}

	// Undoes what remove_sites did, using its returned vector
	void restore_sites(std::vector<bindingsite> &old_sites)
	{
		_sites = old_sites;
	}

	// Check lowest Hamming distance to a TF
	int tfactorhamming(const tfactor &tf) const;

	// Returns the number of 1-bits in cisregion
	int onebits_cis() const
	{
		return _genome.substr(cisregion_pos, cisregion_size).popcount();
	}

	// Returns the number of 1-bits in the binding motif
	int onebits_motif() const
	{
		return _genome.substr(bindmotif_pos, bindmotif_size).popcount();
	}

	// Check _all_ hamming distances to a TF
	void tfactorhammingdists(const tfactor &tf, std::vector<int> &hammings) const;

	// Updates the binding sites using tfactor_count random TFs (that
	// all bind to at least one site).
	void randomize_sites(size_t tfactor_count,
		int maxhamming = maxhamming_dna);

	// Randomizes tfactor_count TFs, saves statistics on the number of binding
	// sites for each of them, and updates the binding sites.
	size_t try_binding(size_t tfactor_count, std::map<int, size_t> &bindSameHist,
		int maxhamming = maxhamming_dna);

	// Return index to all genes in tfactors binding to this cisregion.
	void getbindingtfs(std::set<size_t> &tfindices) const;

	// Return binding start pos to all genes in tfactors binding to this cisregion.
	void getbindingtfspos(std::map<size_t, size_t> &positions) const;

	// Determine whether a TF is activating or repressing transcription. At
	// most one can be true, but neither may be true if it's ambiguous. If
	// the TF has a binding site in the blocking (Buchler: promoter) region,
	// it is considered repressing. Otherwise, if enable_repressors is true
	// we require that at least two thirds of the binding sites agree on the
	// sign for it to be well-defined.
	// Returns positive for activation, negative for repression and 0 for
	// ambiguous/neither.
	int interaction_sign(size_t tfindex) const;

	// Finds the fraction of time that RNAP is bound, using the list
	// of binding sites and the TF levels.
	double getrnapactivity(const std::vector<double> &tflevels) const;

	// An unoptimized implementation of getrnapactivity, purely for
	// testing purposes.
	double getrnapactivity_safe(const std::vector<double> &tflevels) const;

	//save state to out-stream.
	void save(std::ostream &out) const;

	void load(std::istream &in);

	void printbindingsites(std::ostream &out) const;
	void printrnapactivities(std::ostream &out) const;

private:
	static inline bool in_promoter(int pos)
	{
		assert(cisregion_pos <= pos && pos < cisregion_pos + cisregion_size);
		return pos > promoter_start - bindmotif_size && pos < promoter_end;
	}
	inline bool is_repressor(int pos) const
	{
		assert(cisregion_pos <= pos && pos < cisregion_pos + cisregion_size);
		return enable_repressors &&
			!_genome[pos ? pos - 1 : cisregion_pos + cisregion_size - 1];
	}

	// update all parameters from the bitstring representation
	void refresh();

	// Finds and adds binding sites for a specific TF.
	size_t addtfactorsites(const tfactor &tf, size_t index);

	double calcrnapactivity(int tfs) const;

	// The list of binding sites, valid after updatesites()
	std::vector<bindingsite> _sites;

	// RNAP activity as function the of number of bound TFs, tabulated
	// for efficiency.
	std::vector<double> _rnapactivities;

	// Current concentration level of the gene product (protein)
	double _level;

	// Is the protein a tfactor or a protein degrader/protector?
	bool _isregulator;

	// A site on the protein where other proteins can interact, and mark it for degradation.
	struct interaction
	{
		interaction() {}
		size_t proteinindex;	// Index of the protein in the vector of proteins
		double affinity;		// Binding affinity
		bool isdegrader;		// Binds to degrade or protect the protein?
	};

	// The list of interactions, valid after updateinteractions()
	std::vector<interaction> _interactions;

	// The energy by which RNAP binds to DNA (in kT).
	double _rnap_dna_energy;

	// Production rate, to be multiplied with RNAP activity
	double _productionrate;

	// For light (de)activated proteins, the light sensitivity
	double _lightsensitivity;

	// Is the protein activated by light(1) or darkness(-1) or neither(0)?
	int _lightactivation;

	// When active, is it binding to protein to degrade it, or to protect it
	bool _isdegrader;

	// Our full genome bitstring
	bitstring _genome;


public:

	// The binding affinity with DNA as a function of Hamming distance
	static inline double dnaaffinity(int distance)
	{
		// alpha = 1000, affinity decreasing by factor sqrt(10) per mismatch
		return 1000. * std::pow(sqrt(10.), -(double)distance);
	}

	// The binding affinity with proteins as a function of Hamming distance
	static inline double proteinaffinity(int distance)
	{
		// alpha = 1, affinity decreasing by factor sqrt(10) per mismatch
		return 100 * std::pow(sqrt(10.), -(double)distance);
	}

	inline bitstring getbindingmotif() const
	{
		return _genome.substr(bindmotif_pos, bindmotif_size);
	}

	inline double getproductionrate() const
	{
		return _productionrate;
	}

	inline double getdecayrate(const std::vector<double> &proteinlevels) const
	{
		const double base_rate = .01, max_rate = 10;
		double nominator = 0, denominator = 0;
		for(auto & ion : _interactions)
		{
			double tmp = proteinlevels[ion.proteinindex] * ion.affinity;
			denominator += tmp;
			if(ion.isdegrader)
				nominator += tmp;
		}
		return (base_rate + max_rate * nominator) / (1 + denominator);
	}

	inline double getactivity(double light) const
	{
		if(!_lightactivation)
			return 1;
		double a = 1 / (1 + _lightsensitivity * light);
		if(_lightactivation < 0)	// light inactivation?
			return a;
		return 1. - a;
	}

	// Get/set concentration level, to store beteen simulation runs
	inline double getlevel() const
	{
		return _level;
	}
	inline void setlevel(double level)
	{
		_level = level;
	}

	inline bool isregulator() const
	{
		return _isregulator;
	}

	inline bool isdegrader() const
	{
		return _isdegrader;
	}

	inline int get_lightactivation() const
	{
		return _lightactivation;
	}

	inline double get_lightsensitivity() const
	{
		return _lightsensitivity;
	}

	// Return index of all proteins binding to this protein
	void getbindingproteins(std::set<size_t> &binderindices) const
	{
		for(auto &in : _interactions)
			binderindices.insert(in.proteinindex);
	}

	void getbindingproteins(std::set<size_t> &binderindices, bool isdegrader) const
	{
		for(auto &in : _interactions)
			if(isdegrader && in.isdegrader)
				binderindices.insert(in.proteinindex);
			else if(!isdegrader && !in.isdegrader)
				binderindices.insert(in.proteinindex);
	}

	// Removes the interaction, and returns the old vector of interactions
	std::vector<interaction> remove_interaction(size_t index)
	{
		std::vector<interaction> old_interactions(_interactions);
		_interactions.erase(std::remove_if(_interactions.begin(), _interactions.end(),
			[index](interaction x){return index == x.proteinindex;}), _interactions.end());

		return old_interactions;
	}

	// Undoes what remove_sites did, using its returned vector
	void restore_interactions(std::vector<interaction> &old_interactions)
	{
		_interactions = old_interactions;
	}
   
   friend bool operator==(const gene &gene1, const gene &gene2);
   friend bool operator<(const gene &gene1, const gene &gene2);
};


#endif
