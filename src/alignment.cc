#include <algorithm>
#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include <ctime>
#include <exception>
#include "alignment.h"

//#define DEBUGME

static constexpr int gap = -20;
static constexpr int extension = -3;
static constexpr int mismatch = -5;
static constexpr int match = 1;

static inline constexpr int gap_cost(int width)
{
	return gap + (width-1) * extension;
}


#ifdef DEBUGME
static std::ostream &operator<<(std::ostream &out, const std::vector<char> &vec)
{
	for(auto v : vec)
		out << (char)(v < 10 ? v + '0' : v);
	return out;
}

static std::ostream &operator<<(std::ostream &out, const std::vector<int> &vec)
{
	for(auto v : vec)
		out << " " << v;
	return out;
}

static int depth = 0;
std::string spaces()
{
	return std::string(depth, ' ');
}
#define DEB(x) std::cout<<spaces()<<x;
#else
#define DEB(x)
#endif

// Function for aligning according to two scoring vectors in the Needleman-Wunsch matrix
static void lastLineAlign(const std::vector<char> &column, const std::vector<char> &row,
	std::vector<int> &prev, std::vector<int> &prevg, bool initialGap)
{
DEB("LLA "<<column<<" "<<row<<"\n")
DEB(" gap: "<<initialGap<<"\n")
	
	size_t minLen = row.size();
	std::vector<int> cur(minLen + 1);
	std::vector<int> curg(minLen + 1);
	prev = cur;
	prevg = cur;

	prevg[0] = initialGap ? gap + extension: 0;
	prev[0] = 0;

	for(size_t i = 1; i <= minLen; i++)
	{
		prev[i] = initialGap ? i * extension : gap_cost(i);
		prevg[i] = initialGap ? i * extension : gap_cost(i);
		//prev[i] = gap_cost(i);
	}

	for(size_t j = 1; j <= column.size(); j++)
	{
DEB(" prev: "<<prev<<"\n")
DEB(" prevg: "<<prevg<<"\n")
		cur[0] = initialGap ? j * extension : gap_cost(j);
		curg[0] = initialGap ? j * extension : gap_cost(j);

		for(size_t i = 1; i <= minLen; i++)
		{
			// optimal option with gap
			curg[i] = std::max(std::max(cur[i-1], prev[i]) + gap,
				std::max(curg[i-1], prevg[i]) + extension);
			// if match
			if(row[i-1] == column[j-1])
				cur[i] = std::max(prev[i-1], prevg[i-1]) + match;
			else
				cur[i] = std::max(prev[i-1], prevg[i-1]) + mismatch;
		}
		prev.swap(cur);
		prevg.swap(curg);
	}
DEB(" prev: "<<prev<<"\n")
DEB(" prevg: "<<prevg<<"\n")
}

// Function for deciding optimal slicing point for Hirschbergs D&C approach
static int partitionY(const std::vector<int> &scoreL, const std::vector<int> &scoreLg,
	const std::vector<int> &scoreR, const std::vector<int> &scoreRg,
	bool &leftGap, bool &rightGap)
{
DEB("PY "<<leftGap<<" "<<rightGap<<"\n")
DEB(" scoreL "<<scoreL<<"\n")
DEB(" scoreLg "<<scoreLg<<"\n")
DEB(" scoreR "<<scoreR<<"\n")
DEB(" scoreRg "<<scoreRg<<"\n")
	int max_index = -1;
	int max_sum = -99999999;

	for(size_t i = 0; i < scoreL.size(); i++)
	{
		size_t j = scoreL.size()-1-i;

		int best = std::max(scoreL[i] + scoreR[j],
			std::max(scoreLg[i] + scoreR[j],
			std::max(scoreL[i] + scoreRg[j],
			scoreLg[i] + scoreRg[j] - gap + extension)));

		if(best > max_sum)
		{
			max_sum = best;
			max_index = i;

			if(scoreL[i] + scoreR[j] == best)
				leftGap = rightGap = false;
			else if(scoreLg[i] + scoreR[j] == best)
			{
				leftGap = true;
				rightGap = false;
			}
			else if(scoreL[i] + scoreRg[j] == best)
			{
				leftGap = false;
				rightGap = true;
			}
			else
				leftGap = rightGap = true;
		}
	}
	assert(max_index >= 0);
DEB(" max_index "<<max_index<<"\n")
	return max_index;
}

// Builds output as vectors of chars.
static void alignBuilder(const std::vector<char> &column, const std::vector<char> &row,
	std::vector<char>& columnf, std::vector<char>& rowf, bool leftGap, bool rightGap)
{
DEB("AB "<<leftGap<<" "<<rightGap<<"\n")
DEB(" column "<<column<<"\n")
DEB(" row    "<<row<<"\n")
	assert(row.size() == 1);
	if(column.size() == 1)
	{
		char c = column.front();
		char r = row.front();
		// cost of aligning c with r
		int m = r == c ? match : mismatch;
		// cost of aligning c- with -r
		int g = (rightGap ? extension : gap) + (leftGap ? extension : gap);

		columnf.push_back(c);
		if(g >= m)
		{
			columnf.push_back('-');
			rowf.push_back('-');
		}
		rowf.push_back(r);
	}
	else
	{
		char r = row.front();

		std::vector<int> sc(3);
		// align r--- with cccc
		sc[0] = (r == column.front() ? match : mismatch) + (rightGap ? extension : gap); //+(n-2)*ext
		// align ---r with cccc
		sc[1] = (r == column.back() ? match : mismatch) + (leftGap ? extension : gap);  //+(n-2)*ext
		int ngaps = (int)rightGap + leftGap;
		// align r---- with -cccc
		sc[2] = 2 * extension + gap + (extension - gap) * ngaps;  //+(n-2)*ext
		size_t matchpos = 0;
		if(column.size() > 2)
		{
			// we assume match > mismatch
			matchpos = std::find(column.begin() + 1, column.end() - 1, r) - column.begin();
			// align -r-- with cccc
			int smid = 2 * gap - extension + (extension - gap) * ngaps;  //+(n-2)*ext
			sc.push_back(smid + (matchpos < column.size() - 1 ? match : mismatch));
		}
DEB(" ngaps "<<ngaps<<" scores "<<sc<<"\n")
		auto best = std::max_element(sc.begin(), sc.end()) - sc.begin();
		if(best == 0)
		{
			rowf.push_back(r);
			rowf.resize(rowf.size() + column.size() - 1, '-');
		}
		else if(best == 1)
		{
			rowf.resize(rowf.size() + column.size() - 1, '-');
			rowf.push_back(r);
		}
		else if(best == 2)
		{
			rowf.push_back(r);
			rowf.resize(rowf.size() + column.size(), '-');
			columnf.push_back('-');
		}
		else // 3, match/mismatch
		{
			rowf.resize(rowf.size() + matchpos, '-');
			rowf.push_back(r);
			rowf.resize(rowf.size() + column.size() - 1 - matchpos, '-');
		}
		columnf.insert(columnf.end(), column.begin(), column.end());
	}
}

static void Hirschberg(const std::vector<char> &column, const std::vector<char> &row,
	std::vector<char>& columnff, std::vector<char>& rowff, bool leftGap, bool rightGap)
{
#ifdef DEBUGME
depth++;
#endif
DEB("Hirshberg "<<leftGap<<" "<<rightGap<<"\n")
DEB(" column "<<column<<"\n")
DEB(" row "<<row<<"\n")
	if(column.size() == 0 || row.size() == 0) // IF ANY = 0
	{
		if(column.size() == 0)
		{
			columnff.resize(columnff.size() + row.size(), '-');
			rowff.insert(rowff.end(), row.begin(), row.end());
		}
		else
		{
			rowff.resize(rowff.size() + column.size(), '-');
			columnff.insert(columnff.end(), column.begin(), column.end());
		}
	}
	else if(row.size() == 1)
	{
		alignBuilder(column, row, columnff, rowff, leftGap, rightGap);
	}
	else if(column.size() == 1)
	{
		alignBuilder(row, column, rowff, columnff, leftGap, rightGap);
	}
	else
	{
DEB(" recursing\n")
		size_t xlen = column.size();
		size_t xmid = xlen / 2;

		std::vector<char> lcolumn(column.begin(), column.begin() + xmid);
		std::vector<char> rcolumn(column.begin() + xmid, column.end());
		std::vector<int> scoreL, scoreLg, scoreR, scoreRg;

		lastLineAlign(lcolumn, row, scoreL, scoreLg, leftGap);

		std::vector<char> rcolumn_reversed = rcolumn;
		std::vector<char> row_reversed = row;
		std::reverse(rcolumn_reversed.begin(), rcolumn_reversed.end()); //flipping rcol
		std::reverse(row_reversed.begin(), row_reversed.end());	//flipping row

		lastLineAlign(rcolumn_reversed, row_reversed, scoreR, scoreRg, rightGap);

		bool subLeftGap = false, subRightGap = false;
		int ymid = partitionY(scoreL, scoreLg, scoreR, scoreRg, subLeftGap, subRightGap);

		std::vector<char> lrow(row.begin(), row.begin() + ymid);
		std::vector<char> rrow(row.begin() + ymid, row.end());
		std::vector<char> row_l, column_u;
		std::vector<char> row_r, column_d;
		Hirschberg(lcolumn, lrow, column_u, row_l, leftGap, subLeftGap);
		Hirschberg(rcolumn, rrow, column_d, row_r, subRightGap, rightGap);

		rowff.insert(rowff.end(), row_l.begin(), row_l.end());
		rowff.insert(rowff.end(), row_r.begin(), row_r.end());
		columnff.insert(columnff.end(), column_u.begin(), column_u.end());
		columnff.insert(columnff.end(), column_d.begin(), column_d.end());
	}
DEB(" hirshberg result\n")
DEB(" columnff "<<columnff<<"\n")
DEB(" rowff    "<<rowff<<"\n")
#ifdef DEBUGME
depth--;
#endif
}

void align_hirschberg_nc(const bitstring &ba, const bitstring &bb,
	std::vector<char> &column, std::vector<char> &row)
{
	std::vector<char> v(ba.bits());
	std::vector<char> w(bb.bits());

	for(size_t i = 0; i < ba.bits(); i++)
	{
		if(ba[i])
			v[i] = 1;
	}
	for(size_t i = 0; i < bb.bits(); i++)
	{
		if(bb[i])
			w[i] = 1;
	}
	Hirschberg(v, w, column, row, false, false);
}

// Clean up ugly things like 11-00- aligned with 111--1
void tidy_alignment(std::vector<char> &column, std::vector<char> &row)
{
	assert(row.size() == column.size());
	std::vector<char> tmp;
	bool first_rgap = false, first_cgap = false;
	for(size_t i = 0; i < row.size(); i++)
	{
		bool cgap = column[i] == '-', rgap = row[i] == '-';
		if(!cgap && !rgap)
		{
			if(!tmp.empty())
			{
				std::copy(tmp.begin(), tmp.end(),
					(first_rgap ? row : column).begin() + i - tmp.size());
				std::fill_n((first_rgap ? column : row).begin() + i - tmp.size(),
					tmp.size(), '-');
				tmp.clear();
			}
			first_rgap = false;
			first_cgap = false;
		}
		else if(first_rgap || first_cgap)
		{
			if(rgap != first_rgap)
				tmp.push_back((rgap ? column : row)[i]);
			else
			{
				row[i - tmp.size()] = row[i];
				column[i - tmp.size()] = column[i];
			}
		}
		else
		{
			first_rgap = rgap;
			first_cgap = cgap;
		}
	}
	if(!tmp.empty())
	{
		std::copy(tmp.begin(), tmp.end(),
			(first_rgap ? row : column).end() - tmp.size());
		std::fill_n((first_rgap ? column : row).end() - tmp.size(),
			tmp.size(), '-');
	}
}

void align_hirschberg(const bitstring &ba, const bitstring &bb,
	std::vector<char> &column, std::vector<char> &row)
{
	column.clear();
	row.clear();
	align_hirschberg_nc(ba, bb, column, row);
	tidy_alignment(column, row);
}

// Divide and conquer algorithm for reducing computation time if both
// strings share global similarities.
int align_heuristic_rec(const bitstring &ba, const bitstring &bb,
	std::vector<char> &finalc, std::vector<char> &finalr, const int maxdiff)
{
	// Length of substring
	const size_t len = bitstring::Tbits;

	// If this short, don't bother.
	const size_t minlen = 256;

	if(ba.bits() < minlen || bb.bits() < minlen)
	{
		align_hirschberg_nc(ba, bb, finalc, finalr);
		return 0;
	}

	std::vector<size_t> breaka, breakb;
	const size_t cuts = 3;

	// Check fit to substrings in the different sections
	for(size_t i = 0; i < cuts; i++)
	{
		size_t cutpoint = ba.bits() * (i + 1) / (cuts + 1) - len / 2;
		bitstring sub = ba.substr(cutpoint, len);
		int bestdiff = maxdiff + 1;
		size_t bestj(-1);

		for(size_t j = 0; j <= bb.bits() - sub.bits(); j++)
		{
			int diff = bb.distance(sub, j);
			if(diff < bestdiff)
			{
				bestdiff = diff;
				bestj = j;
			}
			else if(diff == bestdiff)
				bestj = size_t(-1);
		}
		if(bestj != size_t(-1))
		{
			breaka.push_back(cutpoint + len / 2);
			breakb.push_back(bestj + len / 2);
		}
	}

	// Too few unique matches: brute force align
	if(breaka.size() < 2)
	{
		align_hirschberg_nc(ba, bb, finalc, finalr);
		return 0;
	}

	std::vector<size_t> brka, brkb;
	brka.push_back(0);
	brkb.push_back(0);

	// Matches neatly ordered? Don't allow overlapping matches.
	for(size_t i = 0; i < breakb.size(); i++)
	{
		if(i > 0 && breakb[i - 1] + len > breakb[i])
			continue;
		if(i + 1 < breakb.size() && breakb[i] + len > breakb[i + 1])
			continue;
		brka.push_back(breaka[i]);
		brkb.push_back(breakb[i]);
	}

	if(brka.size() < 2)
	{
		align_hirschberg_nc(ba, bb, finalc, finalr);
		return 0;
	}
	// Add the ends to the lists of breaks
	brka.push_back(ba.bits());
	brkb.push_back(bb.bits());

	int madecuts = 1;
	for(size_t i = 1; i < brka.size(); i++)
	{
		madecuts += align_heuristic_rec(
			ba.substr(brka[i-1], brka[i] - brka[i-1]),
			bb.substr(brkb[i-1], brkb[i] - brkb[i-1]),
			finalc, finalr, maxdiff);
	}
	return madecuts;
}

// Divide and conquer algorithm for reducing computation time if both
// strings share global similarities.
int align_heuristic(const bitstring &ba, const bitstring &bb,
	std::vector<char> &column, std::vector<char> &row, const int maxdiff)
{
	column.clear();
	row.clear();
	auto tmp = align_heuristic_rec(ba, bb, column, row, maxdiff);
	tidy_alignment(column, row);
	return tmp;
}


// Divide and conquer algorithm for reducing computation time if both
// strings share global similarities.
void align_heuristic_new_rec(const bitstring &ba, const bitstring &bb,
	std::vector<char> &finalc, std::vector<char> &finalr, const int maxdiff,
	const int minmargin, const int retries)
{
	// Length of substring
	const size_t len = bitstring::Tbits;

	// If this short, don't bother.
	const size_t minlen = len * 2;

	if(ba.bits() < minlen || bb.bits() < minlen)
	{
		align_hirschberg_nc(ba, bb, finalc, finalr);
		return;
	}

	for(int i = 0; i < retries; i++)
	{
		size_t cutpoint = i ? gsl_rng_uniform_int(rng, ba.bits() - len) :
			(ba.bits() - len) / 2;
		bitstring sub = ba.substr(cutpoint, len);
		int bestdiff = maxdiff + 1;
		int nextbest = bestdiff;
		size_t bestj(-1);

		for(size_t j = 0; j <= bb.bits() - sub.bits(); j++)
		{
			int diff = bb.distance(sub, j);
			if(diff <= bestdiff)
				nextbest = bestdiff;
			if(diff < bestdiff)
			{
				bestdiff = diff;
				bestj = j;
			}
		}

		if(nextbest - bestdiff >= minmargin)
		{
			size_t cp1 = cutpoint + len / 2, cp2 = bestj + len / 2;
			align_heuristic_new_rec(ba.substr(0, cp1),
				bb.substr(0, cp2), finalc, finalr,
				maxdiff, minmargin, retries);
			align_heuristic_new_rec(ba.substr(cp1, ba.bits() - cp1),
				bb.substr(cp2, bb.bits() - cp2), finalc, finalr,
				maxdiff, minmargin, retries);
			return;
		}
	}

	align_hirschberg_nc(ba, bb, finalc, finalr);
}

// Divide and conquer algorithm for reducing computation time if both
// strings share global similarities.
void align_heuristic_new(const bitstring &ba, const bitstring &bb,
	std::vector<char> &column, std::vector<char> &row, const int maxdiff,
	const int minmargin, const int retries)
{
	column.clear();
	row.clear();
	align_heuristic_new_rec(ba, bb, column, row,
		maxdiff, minmargin, retries);
	tidy_alignment(column, row);
}


class substring
{
public:
	substring(const bitstring *str, size_t begin, size_t size)
		: parent(str), start(begin), length(size) {}

	substring() : parent(0), start(0), length(0) {}

	bitstring get() const{
		return parent->substr(start, length);
	}

	bool operator<(const substring &o) const
	{
		return parent->less_substr(start, length, *o.parent, o.start, o.length);
	}

	// Number of identical bits starting from start
	size_t start_similarity(const substring &o) const
	{
		size_t minlen = std::min(length, o.length);
		for(size_t i = 0, j = o.start; i < minlen; ++i, ++j)
		{
			if((*parent)[i + start] != (*o.parent)[j])
				return i;
		}
		return minlen;
	}

	const bitstring *parent;
	size_t start;
	size_t length;
};



// Find the longest common subsequence (LCSS) of two sequences.
void lcss_sorting(const bitstring &a, const bitstring &b,
	substring &asuffix, substring &bsuffix)
{
	// Fill with all possible suffixes from both a, b.
	std::vector<substring> suffixes;

	for(size_t i = 0; i < a.bits(); ++i)
		suffixes.push_back(substring(&a, i, a.bits() - i));

	for(size_t i = 0; i < b.bits(); ++i)
		suffixes.push_back(substring(&b, i, b.bits() - i));

	std::sort(suffixes.begin(), suffixes.end());

	// Now the two substrings starting at the LCSS (longest common
	// subsequence) must be next to each other in the list. Find the
	// two neighboring substrings with the highest start similarity
	int maxsim = 0;

	asuffix = substring(&a, 0, 0);
	bsuffix = substring(&b, 0, 0);

	for(size_t i = 0; i + 1 < suffixes.size(); i++)
	{
		// Don't compare to suffix in the same sequence
		if(suffixes[i].parent == suffixes[i+1].parent)
			continue;

		int sim = suffixes[i].start_similarity(suffixes[i+1]);
		if(sim > maxsim)
		{
			maxsim = sim;
			if(suffixes[i].parent == &a)
			{
				asuffix.start = suffixes[i].start;
				bsuffix.start = suffixes[i+1].start;
			}
			else
			{
				bsuffix.start = suffixes[i].start;
				asuffix.start = suffixes[i+1].start;
			}
		}
	}
	asuffix.length = bsuffix.length = maxsim;
}

static void append(std::vector<char> &vec, const bitstring &b)
{
	for(size_t i = 0; i < b.bits(); ++i)
		vec.push_back(b[i]);
}

void align_synapsing_rec(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row, size_t tresholdlength,
	size_t align_tresholdlength)
{
	substring lcss_a, lcss_b;
	lcss_sorting(a, b, lcss_a, lcss_b);
	if(lcss_a.length < tresholdlength)
	{
		append(column, a);
		if(a.bits() > align_tresholdlength || a.bits() != b.bits())
		{
			column.resize(column.size() + b.bits(), '-');
			row.resize(row.size() + a.bits(), '-');
		}
		append(row, b);
		return;
	}

	bitstring aleft = a.substr(0, lcss_a.start);
	bitstring bleft = b.substr(0, lcss_b.start);
	align_synapsing_rec(aleft, bleft, column, row, tresholdlength, align_tresholdlength);

	bitstring mid = a.substr(lcss_a.start, lcss_a.length);
	append(row, mid);
	append(column, mid);

	bitstring aright = a.substr(lcss_a.start + lcss_a.length,
		a.bits() - (lcss_a.start + lcss_a.length));
	bitstring bright = b.substr(lcss_b.start + lcss_b.length,
		b.bits() - (lcss_b.start + lcss_b.length));
	align_synapsing_rec(aright, bright, column, row, tresholdlength, align_tresholdlength);
}


void align_synapsing(const bitstring &a, const bitstring &b, std::vector<char> &column,
	std::vector<char> &row, size_t tresholdlength, size_t align_tresholdlength)
{
	column.clear();
	row.clear();
	align_synapsing_rec(a, b, column, row, tresholdlength, align_tresholdlength);
}


// Builds a fake alignment for one-point crossover:
// aaaa---aaaaaaa
// bbbbbb---bbbbb
// where the crossover point has - on both strands.
void align_onepoint(const bitstring &a, const bitstring &b,
   std::vector<char> &column, std::vector<char> &row,
   size_t pos1, size_t pos2)
{
	column.clear();
	row.clear();

	size_t rev1 = a.bits() - pos1;
	size_t rev2 = b.bits() - pos2;
	append(column, a.substr(0, pos1));
	append(row, b.substr(0, pos2));
	column.resize(std::max(pos1, pos2) + 1 + std::max(rev1, rev2) - rev1, '-');
	row.resize(std::max(pos1, pos2) + 1 + std::max(rev1, rev2) - rev2, '-');
	append(column, a.substr(pos1, rev1));
	append(row, b.substr(pos2, rev2));
}


void align_cutnsplice(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row)
{
	size_t pos1 = gsl_rng_uniform_int(rng, a.bits());
	size_t pos2 = gsl_rng_uniform_int(rng, b.bits());

	align_onepoint(a, b, column, row, pos1, pos2);
}

void saga_algb(const bitstring &a, const bitstring &b, int as, int ae,
	std::vector<int> &out)
{
	std::vector<int> temp(b.bits() + 1);
	out = temp;
	if(as < ae)
	{
		for(int i = as; i < ae; ++i)
		{
			for(size_t j = 0; j < b.bits(); ++j)
				temp[j + 1] = a[i] == b[j] ? out[j] + 1 : std::max(temp[j], out[j + 1]);
			out.swap(temp);
		}
	}
	else
	{
		for(int i = as; i > ae; --i)
		{
			for(size_t j = 0; j < b.bits(); ++j)
			{
				temp[j + 1] = a[i] == b[b.bits() - 1 - j] ? out[j] + 1 :
					std::max(temp[j], out[j + 1]);
			}
			out.swap(temp);
		}
	}
}

void find_crossover_point_saga(const bitstring &a, const bitstring &b,
	size_t &pos1, size_t &pos2)
{
	pos1 = gsl_rng_uniform_int(rng, a.bits());
	std::vector<int> l1, l2;
	saga_algb(a, b, 0, pos1, l1);
	saga_algb(a, b, a.bits() - 1, pos1 - 1, l2);
	int best = 0;
	std::vector<size_t> bests;
	for(size_t i = 0; i < l1.size(); ++i)
	{
		int score = l1[i] + l2[l2.size() - 1 - i];
		if(score > best)
		{
			bests.clear();
			best = score;
		}
		if(score == best)
			bests.push_back(i);
	}
	pos2 = bests[gsl_rng_uniform_int(rng, bests.size())];
}

void align_saga(const bitstring &a, const bitstring &b,
	std::vector<char> &column, std::vector<char> &row)
{
	size_t pos1, pos2;
	find_crossover_point_saga(a, b, pos1, pos2);
	align_onepoint(a, b, column, row, pos1, pos2);
}



int score_alignment(std::vector<char> &seq1, std::vector<char> &seq2)
{
	assert(seq1.size() == seq2.size());
	int score = 0;
	size_t n = seq1.size();
	bool prevgap = false;
	for(size_t i = 0; i < n; ++i)
	{
		auto c = seq1[i], d = seq2[i];
		if(c == '-' && d == '-')
			continue;
		else if(c == '-' || d == '-')
		{
			score += prevgap ? extension : gap;
			prevgap = true;
		}
		else
		{
			score += c == d ? match : mismatch;
			prevgap = false;
		}
	}
	return score;
}


