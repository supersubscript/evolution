#ifndef DIGORG_BOOLEAN_H
#define DIGORG_BOOLEAN_H

#include <vector>
#include <algorithm>

int count_inputs(const std::vector<int> &rule);

// For a given Boolean rule, finds the number of unused
// inputs, and removes them.
// Input: a truth table with 2**n elements for some n.
void prune_unused_input(const std::vector<int> &allrules,
	std::vector<int> &used_rules, int &inp_unused);

// For a given used Boolean rule, finds the number of nested
// canalizing inputs.
// Input: a truth table with 2**n elements for some n, of used rules.
void check_canalization(const std::vector<int> &rule,
	std::vector<int> &other_rule,
	int &inp_canal, std::vector<std::pair<bool, bool>> &canal_inout);

#endif
