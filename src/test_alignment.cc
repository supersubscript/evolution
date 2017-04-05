#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include <gsl/gsl_randist.h>
#include "alignment.cc"

gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);

void pr(const bitstring &b)
{
	int c = 0;
	std::cout << "c\t" ;
	for(size_t i = 0; i < b.bits(); ++i)
	{
		if(c==10) c = 0;
		std::cout << c;
		++c;
	}
	std::cout << std::endl;
}

void align_and_print(bitstring &a, bitstring &b)
{
	std::vector<char> avec, bvec;
	align_hirschberg(a, b, avec, bvec);
	std::cout << "hirschberg:\n";
	for(auto c : avec)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\n";
	for(char c : bvec)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\nscore: " << score_alignment(avec, bvec) << "\n";

	std::cout << "synapsing:\n";
	align_synapsing(a, b, avec, bvec, 3, 0);
	for(auto c : avec)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\n";
	for(char c : bvec)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\nscore: " << score_alignment(avec, bvec) << "\n";

	std::cout << "saga:\n";
	align_saga(a, b, avec, bvec);
	for(auto c : avec)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\n";
	for(char c : bvec)
		std::cout << (char)(c > 2 ? c : c + '0');
	std::cout << "\nscore: " << score_alignment(avec, bvec) << "\n";
}

bool tiny_gap(const std::vector<char> &vec)
{
	int g = 10;
	for(auto c : vec)
	{
		if(c == '-' && g > 0 && g < 3) return true;
		else if(c == '-') g = 0;
		else ++g;
	}
	return false;
}

void tinytest()
{
//	bitstring a("011101110001111110110110100101111111100100101111010011100000001011111110010100001101010011111101101000001111000011111110101111110011110101010011101011000010111001100111111101111011001010011100110000111010100110010011100000101111000110011010001001011001000100110101111111101101100001100100100011011010010110110011101001010111010010111111001011100010011100111001101010101111000010110011001001011001110101110101100101000101101001111101101010111000110000110"),
//		b("11101100011011101110001000101111110000011111111111001010000111101001111000110010101110111111010000111101010001111101001110111010011111001111110011101011000010111001100111111101111011001010011100110000111010000110010011100000001111000110011010111001011001000100000101110111111101100001100101100011010010010110110011101001010111010110111111001001100000011100111001101010100101000111110011001001011001110101110101100101100101101001111101101010111000010000111");
	bitstring a("000100101100100010010101111110110110110000110010010001101101001011011001110100101011101011100110101101111000010111001100100101100111010111010"),
				 b("000100101100100010011010111111110110110000110010010001101101001011011001110100101011101011100110101010111100001011001100100101100111010111010");
	a.print(std::cout);
	std::cout << "\n";
	b.print(std::cout);
	std::cout << "\n";
	align_and_print(a, b);
}

void dostuff()
{
	gsl_rng_set(rng, randrandseed());

//tinytest();
//return;

	int N = 20;

	std::cout << "test 1:\n";
	bitstring a(N);
	a.rand();
	bitstring b = a;
	for(int i = N * 2 / 5; i < N * 3 / 5; ++i)
	{
//		a.set(i, false);
//		b.set(i, true);
	}
	align_and_print(a, b);

	std::cout << "\ntest 2:\n";
	a.rand();
	b = a;
	b.mutate(.1);
	align_and_print(a, b);

	std::cout << "\ntest 3:\n";
	a.rand();
	b = a;
	b.indel_constant_size(-1.5);
	align_and_print(a, b);

	std::cout << "\ntest gap:\n";
	N = 12;
	// TEST 
	for(int j = 0; j < 10000; ++j)
	{
		bitstring a(N), b(N);
		a.rand();
		b.rand();
		std::vector<char> avec, bvec;
		align_hirschberg(a, b, avec, bvec);
		if(tiny_gap(avec) || tiny_gap(bvec))
		{
			std::cout << "\nfound tiny gap:\n";
			align_and_print(a, b);
			break;
		}
	}

	std::cout << "\ntest 4:\n";
	N = 1000;
	// TEST 
	for(int j = 0; j < 1000; ++j)
	{
		bitstring a(N), b(N);
		a.rand();
		b.rand();
		substring lcss_a, lcss_b;
		lcss_sorting(a, b, lcss_a, lcss_b);
	}
}

int main()
{
	try
	{
		dostuff();
	}
	catch(std::string s)
	{
		std::cerr << "error: " << s << "\n";
		return 1;
	}
}
