#include "common.h"
#include "boolean.h"


// Returns true if the rule is canalizing on the given input.
// If true, canal_input and canal_output are set to the canalizing
// input and output. The rule is assumed to depend on all inputs.
static bool is_canalizing(const std::vector<int> &rule, int input,
	bool &canal_input, bool &canal_output)
{
	// The number of outputs that are true when the input is true/false
	size_t on_ons = 0, off_ons = 0;
	size_t bit = size_t(1) << input;

	for(size_t j = 0; j < rule.size(); ++j)
	{
		if(rule[j])
		{
			if(j & bit)
				++on_ons;
			else
				++off_ons;
		}
	}

	if(on_ons == 0)
	{
		canal_input = true;
		canal_output = false;
	}
	else if(on_ons == rule.size() / 2)
	{
		canal_input = true;
		canal_output = true;
	}
	else if(off_ons == 0)
	{
		canal_input = false;
		canal_output = false;
	}
	else if(off_ons == rule.size() / 2)
	{
		canal_input = false;
		canal_output = true;
	}
	else
		return false;
	// else-else
	return true;
}

// Returns the depth to which this rule is canalizing. The
// rule is updated to contain only the non-canalizing part.
static int canalizing_depth(std::vector<int> &rule,
	std::vector<std::pair<bool, bool>> &canal_inout)
{
	if(rule.size() < 2)
		return 0;

	int n = onebits(rule.size() - 1);
	for(int i = 0; i < n; ++i)
	{
		bool canal_input, canal_output;
		if(is_canalizing(rule, i, canal_input, canal_output))
		{
			size_t bit = size_t(1) << i;
			// Construct other_rule from the rule with this input
			// set to its non-canalizing value.
			std::vector<int> newrule;
			for(size_t j = 0; j < rule.size(); ++j)
			{
				if((bool)(j & bit) != canal_input)
					newrule.push_back(rule[j]);
			}
			rule.swap(newrule);
			canal_inout.push_back(std::make_pair(canal_input, canal_output));
			return 1 + canalizing_depth(rule, canal_inout);
		}
	}
	return 0;
}

int count_inputs(const std::vector<int> &rule)
{
	// Find the number of inputs
	int n = 0;
	while(size_t(1) << n < rule.size())
		++n;
	if(size_t(1) << n != rule.size())
		throw std::string("input to check_canalization must have 2**n elements");
	return n;
}

// For a given Boolean rule, finds the number of unused
// inputs, and removes them.
// Input: a truth table with 2**n elements for some n.
void prune_unused_input(const std::vector<int> &rule_all_inp,
	std::vector<int> &rule_used_inp, int &inp_unused)
{
	int n = count_inputs(rule_all_inp);

	// Special case: no inputs.
	if(!n)
	{
		rule_used_inp = rule_all_inp;
		inp_unused = 0;
		return;
	}

	// Unused inputs, stored so we can remove them from the rule
	std::vector<int> unused_inputs;

	// For each input, check if it is unused
	for(int i = 0; i < n; ++i)
	{
		// Is the input used at all?
		bool isUsed = false;

		// Compare all pairs of indexes that differ only by this input's bit.
		size_t bit = size_t(1) << i;

		for(size_t j = 0; j < rule_all_inp.size(); ++j)
		{
			// Only compare each pair once.
			if(j & bit)
				continue;
			// Values with this input off/on
			bool voff = rule_all_inp[j], von = rule_all_inp[j | bit];
			if(voff != von)
				isUsed = true;
		}
		if(!isUsed)
			unused_inputs.push_back(i);
	}
	inp_unused = unused_inputs.size();

	// Describe the rule in terms of used inputs only
	rule_used_inp.clear();
	for(size_t j = 0; j < rule_all_inp.size(); ++j)
	{
		// Extract from truth table only where all unused inputs are false
		bool use = true;
		for(int i : unused_inputs)
		{
			if(j & (size_t(1) << i))
				use = false;
		}
		if(use)
			rule_used_inp.push_back(rule_all_inp[j]);
	}
}

// For a given used Boolean rule, finds the number of nested
// canalizing inputs.
// Input: a truth table with 2**n elements for some n, of used rules.
void check_canalization(const std::vector<int> &rule,
	std::vector<int> &other_rule,
	int &inp_canal, std::vector<std::pair<bool, bool>> &canal_inout)
{
	canal_inout.clear();

	// Find, count and remove canalizing inputs
	other_rule = rule;
	inp_canal = canalizing_depth(other_rule, canal_inout);
}
