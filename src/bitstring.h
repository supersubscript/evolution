#ifndef META_BITSTRING_H
#define META_BITSTRING_H

#include <cassert>
#include <limits>
#include <cinttypes>
#include <vector>
#include <string>
#include <gsl/gsl_randist.h>
#include "common.h"

/*
	A string of 64*n bits
*/

struct bitstring
{
	typedef uint64_t T;

	static constexpr size_t Tbits = std::numeric_limits<T>::digits;

	bitstring()
	{
	}

	// We make _data a bit bigger than needed to have some 0-padding that
	// makes operations easier.
	explicit bitstring(size_t bits) : _bits(bits), _data(usedTs() + 1)
	{
	}

	explicit bitstring(const std::string &s): bitstring(s.size())
	{
		for(size_t i = 0; i < s.size(); ++i)
		{
			if(s[i] == '1')
				set(i, 1);
			else if(s[i] != '0')
				throw "bitstring binary format error in " + s;
		}
	}

	explicit bitstring(const std::vector<char> &v): bitstring(v.size())
	{
		for(size_t i = 0; i < v.size(); ++i)
		{
			if(v[i] == 1)
				set(i, 1);
			else if(v[i] != 0)
				throw std::string("bitstring vector format error ") + v[i];
		}
	}

	void reinitialize(size_t bits)
	{
		_bits = bits;
		_data.resize(usedTs() + 1);
	}

	// Number of bits represented by the bitstring
	inline size_t bits() const
	{
		return _bits;
	}

	inline size_t usedTs() const
	{
		return (_bits + Tbits - 1) / Tbits;
	}

	static inline constexpr T lowmask(int bits)
	{
		return T(-1) >> ((-bits) % Tbits);
	}

	// Randomize entire bitstring
	inline void rand()
	{
		for(size_t i = 0; i < usedTs(); i++)
			_data[i] = T(gsl_rng_get(rng)) << 32 | T(gsl_rng_get(rng));

		// zero out unused trailing bits
		_data[_bits / Tbits] &= lowmask(_bits);
	}

	// Mutates a fraction of the bits
	void mutate(double mutrate)
	{
/*		// Flip bits
		for(size_t p = gsl_ran_geometric(rng, mutrate) - 1; p < _bits;
			p += gsl_ran_geometric(rng, mutrate))
		{
			toggle(p);
		}*/
		int muts = (int)(mutrate * _bits);
		if(gsl_rng_uniform(rng) < mutrate * _bits - muts)
			++muts;
		for(int i = 0; i < muts; ++i)
			toggle(gsl_rng_uniform_int(rng, _bits));
	}

	// Inserts a duplication of a random part of the bitstring, with length
	// drawn from a power law distribution.
	void duplication(double length_power, size_t maxlen)
	{
		assert(_bits < maxlen);
		size_t len = ran_discrete_power(rng, length_power,
			std::min(_bits, maxlen - _bits));
		size_t from = gsl_rng_uniform_int(rng, _bits - len + 1);
		size_t to = gsl_rng_uniform_int(rng, _bits + 1);
		insert(to, substr(from, len));
	}

	// Deletes a random part of the bitstring, with length
	// drawn from a power law distribution.
	void deletion(double length_power, size_t minlen)
	{
		if(_bits <= minlen)
			return;
		size_t len = ran_discrete_power(rng, length_power, _bits - minlen);
		size_t start = gsl_rng_uniform_int(rng, _bits - len);
		remove(start, len);
	}

	// Duplicates a part of the bitstring to a random position and removes
	// the same number of bits at random position.
	// The length is drawn from a discrete power law distribution with
	// indel_length_power <= 0.
	void indel_constant_size(double indel_length_power)
	{
		// duplicate/insert AND delete
		size_t len = ran_discrete_power(rng, indel_length_power, _bits);
		size_t from = gsl_rng_uniform_int(rng, _bits - len + 1);
		size_t to = gsl_rng_uniform_int(rng, _bits + 1);
		insert(to, substr(from, len));
		size_t start = gsl_rng_uniform_int(rng, _bits - len);
		remove(start, len);
	}

	// Returns the number of 1-bits
	inline size_t popcount() const
	{
		size_t p = 0;
		for(size_t i = 0; i < usedTs(); ++i)
			p += onebits(_data[i]);
		return p;
	}

	// Set one bit to 0 or 1
	inline void set(size_t bit, bool val)
	{
		assert(bit < _bits);
		size_t B = bit / Tbits;
		int b = bit % Tbits;
		if(val)
			_data[B] |= T(1) << b;
		else
			_data[B] &= ~(T(1) << b);
	}
	// Toggle one bit
	inline void toggle(size_t bit)
	{
		assert(bit < _bits);
		size_t B = bit / Tbits;
		int b = bit % Tbits;
		_data[B] ^= T(1) << b;
	}

	// Get value of the bit at index 'bit'
	bool operator[](size_t bit) const
	{
		assert(bit < _bits);
		size_t B = bit / Tbits;
		int b = bit % Tbits;
		return _data[B] & T(1) << b;
	}

	bool operator==(const bitstring &o) const
	{
		if(_bits != o.bits())
			return false;
		if(distance(o, 0, 0) == -1)
			return false;
		return true;
	}

	// Get a substring of bits as a number
	inline T get(size_t startpos, size_t len) const
	{
		assert(startpos + len <= _bits);
		assert(len <= Tbits);

		int offs = startpos % Tbits;
		size_t w = startpos / Tbits;
		T ret = _data[w] >> offs;
		if(offs)
			ret |= _data[w + 1] << (Tbits - offs);
		ret &= lowmask(len);
		return ret;
	}

	inline bitstring substr(size_t startpos, size_t len) const
	{
		assert(startpos + len <= _bits);
		bitstring ret(len);
		if(len == 0)
			return ret;

		int offs = startpos % Tbits;
		int w = startpos / Tbits;

		if(!offs)
		{
			for(size_t i = 0; i < ret.usedTs(); ++i, ++w)
				ret._data[i] = _data[w];
		}
		else
		{
			for(size_t i = 0; i < ret.usedTs(); ++i, ++w)
				ret._data[i] = _data[w] >> offs | _data[w + 1] << (Tbits - offs);
		}
		// zero out unused trailing bits
		ret._data[len / Tbits] &= lowmask(len);
		return ret;
	}

	// Add bitstring to end of this
	void append(const bitstring &bitstr)
	{
		size_t w = _bits / Tbits;
		size_t offs = _bits % Tbits;
		_bits += bitstr._bits;
		_data.resize(usedTs() + 1);

		if(!offs)
		{
			for(size_t i = 0; i < bitstr.usedTs(); ++i, ++w)
				_data[w] = bitstr._data[i];
		}
		else
		{
			for(size_t i = 0; i < bitstr.usedTs(); ++i, ++w)
			{
				_data[w] |= bitstr._data[i] << offs;
				_data[w + 1] |= bitstr._data[i] >> (Tbits - offs);
			}
		}
	}

	void truncate(size_t len)
	{
		assert(len <= _bits);
		_bits = len;
		_data.resize(usedTs());
		_data.back() &= lowmask(_bits);
		_data.push_back(T(0));
	}

	void insert(size_t startpos, const bitstring &bitstr)
	{
		assert(startpos <= _bits);
		bitstring b(substr(startpos, _bits - startpos));
		truncate(startpos);
		append(bitstr);
		append(b);
	}


	void remove(size_t startpos, size_t len)
	{
		assert(startpos + len <= _bits);
		bitstring b(substr(startpos + len, _bits - startpos - len));
		truncate(startpos);
		append(b);
	}


	// Returns the Hamming distance to bitstring m when its starting
	// bit is aligned with our bit number pos.
	// pos must be between 0 and bits-m.bits
	// Returns int; broken for strings longer than 2^31-1
	int distance(const bitstring &m, size_t pos) const
	{
		assert(m.bits() + pos <= _bits);
		assert(m.bits() > 0);

		// Split position into words and bits
		int sh = pos % Tbits;
		size_t w = pos / Tbits;
		int d = 0;

		// Are the Ts neatly aligned?
		if(!sh)
		{
			for(size_t i = 0; i < m.usedTs() - 1; ++i, ++w)
				d += onebits(_data[w] ^ m._data[i]);
			d += onebits((_data[w] ^ m._data[m.usedTs() - 1]) & lowmask(m._bits));
		}
		else
		{
			// General case, not neatly aligned.
			for(size_t i = 0; i < m.usedTs() - 1; ++i, ++w)
			{
				d += onebits(m._data[i] ^
					(_data[w] >> sh | _data[w + 1] << (Tbits - sh)));
			}
			T lastdata = _data[w] >> sh | _data[w + 1] << (Tbits - sh);
			d += onebits((m._data[m.usedTs() - 1] ^ lastdata) & lowmask(m._bits));
		}
		return d;
	}
	
   int distance(const bitstring &m, size_t pos, int max_hamming) const
	{
		assert(m.bits() + pos <= _bits);
		assert(m.bits() > 0);

		// Split position into words and bits
		int sh = pos % Tbits;
		size_t w = pos / Tbits;
		int d = 0;

		// Are the Ts neatly aligned?
		if(!sh)
		{
			for(size_t i = 0; i < m.usedTs() - 1; ++i, ++w)
			{	
            d += onebits(_data[w] ^ m._data[i]);
			   if (d > max_hamming)
               return -1;
         }
         d += onebits((_data[w] ^ m._data[m.usedTs() - 1]) & lowmask(m._bits));
	      if (d > max_hamming)
            return -1;
      
      }
		else
		{
			// General case, not neatly aligned.
			for(size_t i = 0; i < m.usedTs() - 1; ++i, ++w)
			{
				d += onebits(m._data[i] ^
					(_data[w] >> sh | _data[w + 1] << (Tbits - sh)));
			   if (d > max_hamming)
               return -1;
         }
			T lastdata = _data[w] >> sh | _data[w + 1] << (Tbits - sh);
			d += onebits((m._data[m.usedTs() - 1] ^ lastdata) & lowmask(m._bits));
         
         if (d > max_hamming)
            return -1;
      }
		return d;
	}

	// Returns the Hamming distance to m when its starting
	// bit is aligned with our bit number pos.
	// pos must be between 0 and bits - Tbits
	int distance(T m, size_t pos) const
	{
		// Split position into words and bits
		int sh = pos % Tbits;
		size_t w = pos / Tbits;

		// Are the Ts neatly aligned?
		if(!sh)
			return onebits(_data[w] ^ m);
		else
			return onebits(m ^ (_data[w] >> sh | _data[w + 1] << (Tbits - sh)));
	}

	// Lexical comparison between a substring of this bitstring and a substring
	// of another. Returns true if this substring is lexically less than the
	// other substring.
	bool less_substr_slow(size_t start, size_t len, const bitstring &o,
		size_t ostart, size_t olen) const
	{
		while(len && olen)
		{
			bool b = o[ostart++];
			if((*this)[start++] != b)
				return b;
			--len;
			--olen;
		}
		return (bool)olen;
	}


	// Lexical comparison between a substring of this bitstring and a substring
	// of another. Returns true if this substring is lexically less than the
	// other substring.
	bool less_substr(size_t start, int len, const bitstring &o,
		size_t ostart, int olen) const
	{
		int offs = start % Tbits, w = start / Tbits;
		int ooffs = ostart % Tbits, ow = ostart / Tbits;

		for(int minlen = std::min(len, olen); minlen > 0; minlen -= Tbits, ++w, ++ow)
		{
			T a = _data[w], b = o._data[ow];
			if(offs)
				a = a >> offs | _data[w + 1] << (Tbits - offs);
			if(ooffs)
				b = b >> ooffs | o._data[ow + 1] << (Tbits - ooffs);
			a ^= b;
			if(a)
			{
				int pos = __builtin_ffsl(a);
				if(pos > minlen)
					return len < olen;
				return b & T(1) << (pos - 1);
			}
		}
		return len < olen;
	}

	/*
	static inline uint32_t reverse(uint32_t v)
	{
		T r = (v & 0x55555555) << 1 | (v & 0xaaaaaaaa) >> 1;
		r = (r & 0x33333333) << 2 | (r & 0xcccccccc) >> 2;
		r = (r & 0x0f0f0f0f) << 4 | (r & 0xf0f0f0f0) >> 4;
		r = (r & 0x00ff00ff) << 8 | (r & 0xff00ff00) >> 8;
		return (r & 0x0000ffff) << 16 | (r & 0xffff0000) >> 16;
	}

	void reverse()
	{
		std::vector<T> rev(usedTs() + 1);
		int sh = _bits % Tbits;
		if(!sh)
		{
			for(size_t i = 0; i < usedTs(); ++i)
				rev[i] = reverse(_data[usedTs() - 1 - i]);
			rev.back() = 0;
		}
		else
		{
			assert(false);
		}
	}*/

	void print(std::ostream &o) const
	{
		for(size_t i = 0; i < _bits; i++)
			o << ((*this)[i] ? "1" : "0");
	}

	void write(std::ostream &out) const
	{
		out << _bits << " ";
		for(size_t b = 0; b < _bits; b += 6)
		{
			int v = 0;
			for(size_t a = 0; a < 6 && a + b < _bits; ++a)
			{
				if((*this)[a + b])
					v |= 1 << a;
			}
			out << (char)('0' + v);
		}
		out << "\n";
	}

	void read(std::istream &in)
	{
		in >> _bits;
		if(!in)
			throw std::string("error reading bit count for bitstring");
		_data.assign(usedTs() + 1, 0);
		if(!_bits)
			return;
		std::string s;
		if(!(in >> s))
			throw std::string("error reading data for bitstring");
		if(s.size() != (_bits + 5) / 6)
			throw std::string("truncated data for bitstring");
		for(size_t b = 0; b < _bits; b += 6)
		{
			int v = s[b / 6] - '0';
			if(v < 0 || v >= 64)
				std::string("invalid encoded bit data for bitstring");
			for(size_t a = 0; a < 6 && a + b < _bits; ++a)
				set(a + b, v & (1 << a));
		}
	}

private:
	size_t _bits;
	std::vector<T> _data;
};

static inline std::ostream &operator<<(std::ostream &o, const bitstring &b)
{
	b.write(o);
	return o;
}

static inline std::istream &operator>>(std::istream &i, bitstring &b)
{
	b.read(i);
	return i;
}


#endif
