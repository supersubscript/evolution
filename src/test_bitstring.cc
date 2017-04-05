#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include "bitstring.h"

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


void dostuff()
{
	gsl_rng_set(rng, randrandseed());

	int N = 100;

	// test << >>
	for(int i = 1; i < N; i++)
	{	
		for(int j = 0; j < 100; ++j)
		{
			bitstring bs(i), bs2;
			bs.rand();
			std::ostringstream s1;
			s1 << bs;
			std::istringstream s2(s1.str());
			s2 >> bs2;
			std::ostringstream s3, s4;
			bs.print(s3);
			bs2.print(s4);
			assert(s3.str() == s4.str());
		}
	}
	std::cout << "<< >> passed\n";

	// test substr
	for(int i = 1; i < N; i++)
	{	
		bitstring bs(i);
		bs.rand();
		std::ostringstream s1;
		bs.print(s1);

		for(int j = 0; j < i; j++)
		{
			for(int k = 1; k < i - j; k++)
			{
				std::ostringstream s2;
				bs.substr(j, k).print(s2);
				assert(s2.str() == s1.str().substr(j, k));
			}
		}
	}
	std::cout << "substr passed\n";


	// TEST append() function
	for(int k = 0; k < N; ++k){
		for(int j = 0; j < N; ++j){
			bitstring b(k);
			bitstring a(j);
			for(size_t i = 0; i < a.bits(); ++i)
				a.set(i, 1);

			std::ostringstream s1,s2,s3;
			b.print(s1);
			a.print(s2);
			b.append(a);
			b.print(s3);

			assert(s1.str()+s2.str() == s3.str());
		}
	}
	std::cout << "append passed\n";

	// TEST less_substr
	N = 200;
	for(int j = 0; j < 10000; ++j){
		bitstring a(1 + gsl_rng_uniform_int(rng, N));
		bitstring b(1 + gsl_rng_uniform_int(rng, N));
		a.rand();
		b.rand();
		for(int k = 0; k < 1000; ++k){
			size_t starta = gsl_rng_uniform_int(rng, a.bits());
			size_t startb = gsl_rng_uniform_int(rng, b.bits());
			size_t lena = gsl_rng_uniform_int(rng, a.bits() - starta + 1);
			size_t lenb = gsl_rng_uniform_int(rng, b.bits() - startb + 1);
			if(a.less_substr(starta, lena, b, startb, lenb) !=
				a.less_substr_slow(starta, lena, b, startb, lenb))
			{
				a.print(std::cerr);
				std::cerr << "\n";
				b.print(std::cerr);
				std::cerr << "\n" << starta << ":" << lena << " vs " <<
					startb << ":" << lenb << "\n";
				a.substr(starta, lena).print(std::cerr);
				std::cerr << "\n";
				b.substr(startb, lenb).print(std::cerr);
				std::cerr << "\n";
				std::cerr << "\n" << a.less_substr(starta, lena, b, startb, lenb) <<
					" != " << a.less_substr_slow(starta, lena, b, startb, lenb) << "\n";
				throw std::string("less_substr failed");
			}
		}
	}
	std::cout << "less_substr_opt passed\n";

	// TEST distance() function
	for(int j = 0; j < 100000; ++j){
		bitstring a(1 + gsl_rng_uniform_int(rng, N));        // substrate
		bitstring b(1 + gsl_rng_uniform_int(rng, a.bits())); // "enzyme"
		a.rand();
		b.rand();

		for(int pos = 0; pos < (int)a.bits() - (int)b.bits() +1; ++pos)
		{
			std::ostringstream s1,s2;
			a.print(s1);
			b.print(s2);

			std::string s_a(s1.str(), pos, s2.str().size());
			std::string s_b(s2.str());

			int sum = 0;
			for(size_t i = 0; i < s_a.size(); ++i)
				if(s_a[i] != s_b[i])
					++sum;
			int x = a.distance(b, pos);
			assert(x==sum);
		}
	}
	std::cout << "distance passed\n";
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
