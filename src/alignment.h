#ifndef DIGORG_ALIGNMENT_H
#define DIGORG_ALIGNMENT_H
#include "bitstring.h"

// Global alignment using the recursive Hirschberg algorithm.
// column and row will contain the aligned a and b
// represented as 0, 1 and '-'. (Note: not '0' and '1'.)
void align_hirschberg(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row);

// Alignment using the "three substrings" heuristic.
// maxdiff is the maximum number of mismatching bits for a substring of
// one sequence to be considered to match some region of the other.
// Returns the total number of cuts (as an indicator of the success of
// the heuristic).
int align_heuristic(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row, const int maxdiff);

// One substring heuristic
void align_heuristic_new(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row, const int maxdiff,
	const int minmargin, const int retries);

// Synapsing alignment.
void align_synapsing(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row, size_t tresholdlength,
	size_t align_tresholdlength);

void align_cutnsplice(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row);

void find_crossover_point_saga(const bitstring &a, const bitstring &b,
	size_t &pos1, size_t &pos2);
void align_saga(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row);

// Scores an alignment according to the costs used by align_hirschberg.
int score_alignment(std::vector<char> &seq1, std::vector<char> &seq2);

#endif
