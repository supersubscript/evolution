#ifndef META_COMMON_H
#define META_COMMON_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <map>
#include <sstream>
#include <cassert>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <cmath>

extern gsl_rng *rng;


static inline double sqr(double v)
{
	return v * v;
}

// returns the number of bits that are set
static inline int onebits(const uint32_t &v)
{
	return __builtin_popcount(v);
}
static inline int onebits(const uint64_t &v)
{
	return __builtin_popcountll(v);
}

// CPU usage since last call, in seconds
static inline double getcputime()
{
	static double cputime = 0;
	rusage rus;
	getrusage(RUSAGE_SELF, &rus);
	double t = double(rus.ru_utime.tv_sec) + double(rus.ru_utime.tv_usec) * 1e-6;
	double diff = t - cputime;
	cputime = t;
	return diff;
}

static inline uint32_t randrandseed()
{
	std::ifstream r("/dev/urandom");
	uint32_t s;
	do
	{
		r.read((char *)&s, sizeof(uint32_t));
	}
	while(!s);
	return s;
}

// p(x) \propto pow(x, power)
static inline size_t ran_discrete_power(gsl_rng *rng, double power, size_t limit)
{
	assert(limit > 0);
	double s = 0;
	for(size_t i = 1; i <= limit; i++)
		s += std::pow(double(i), power);
	s *= gsl_rng_uniform(rng);
	for(size_t i = 1; i <= limit; i++)
	{
		if((s -= std::pow(double(i), power)) <= 0)
			return i;
	}
	return limit;
}

// send vector of zeros, will return all binary permutations or false
// when done.
static inline bool nextbinary(std::vector<int> &onoff)
{
	for(size_t i = 0; i < onoff.size(); ++i)
	{
		if(!onoff[i])
		{
			onoff[i] = 1;
			return true;
		}
		onoff[i] = 0;
	}
	return false;
}

static inline double todouble(const std::string &s)
{
	std::istringstream ss(s);
	double d;
	ss >> d;
	if(!ss || !ss.eof())
		throw std::string("error converting '" + s + "' to double");
	return d;
}

static inline unsigned touint(const std::string &s)
{
	std::istringstream ss(s);
	unsigned u;
	ss >> u;
	if(!ss || !ss.eof())
		throw std::string("'" + s + "' is not a valid unsigned int");
	return u;
}

static inline int toint(const std::string &s)
{
	std::istringstream ss(s);
	int u;
	ss >> u;
	if(!ss || !ss.eof())
		throw std::string("'" + s + "' is not a valid int");
	return u;
}

static inline size_t tosize_t(const std::string &s)
{
	std::istringstream ss(s);
	size_t u;
	ss >> u;
	if(!ss || !ss.eof())
		throw std::string("'" + s + "' is not a valid size_t");
	return u;
}

const std::string &sourcedir(const char *arg = 0);

#endif
