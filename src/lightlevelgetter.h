#ifndef DIGORG_LIGHTLEVELGETTER_H
#define DIGORG_LIGHTLEVELGETTER_H

#include <string>
#include <iostream>
#include <fstream>


// Class that, given a time, returns the light level.
class lightlevelgetter
{
protected:
	lightlevelgetter(double period):
		_period(period)
	{
	}

public:
	virtual ~lightlevelgetter() {}

	double get_period() const
	{
		return _period;
	}

	virtual size_t get_days() const = 0;

	// t is supposed to be between 0 and period
	virtual double get_light(size_t day, double t) const = 0;

	virtual std::string get_name() const = 0;

protected:
	double _period;
};

class lightlevelgetter_ld: public lightlevelgetter
{
public:
	lightlevelgetter_ld(double period, double lightson, double lightsoff):
		lightlevelgetter(period), _lightson(lightson), _lightsoff(lightsoff)
	{
	}

	lightlevelgetter_ld(double period, double photoperiod):
		lightlevelgetter(period),
		_lightson(.5 * (period - photoperiod)),
		_lightsoff(.5 * (period + photoperiod))
	{
	}

	virtual size_t get_days() const
	{
		return 1;
	}

	// t is supposed to be between 0 and period
	virtual double get_light(size_t __attribute__ ((unused)) day, double t) const
	{
		if(t >= _lightson && t <= _lightsoff)
			return 1000;  // because harvard forest data is about this much.
		return 0;
	}

	virtual std::string get_name() const
	{
		std::ostringstream oss;
		double pp = _lightsoff - _lightson;
		oss << "ld" << pp << "_" << _period - pp;
		return oss.str();
	}

private:
	//light can shift twice each period
	double _lightson;		//light on
	double _lightsoff;	//light off
};

class lightlevelgetter_ll: public lightlevelgetter
{
public:
	lightlevelgetter_ll(double period, double level):
		lightlevelgetter(period), _level(level)
	{
	}

	virtual size_t get_days() const
	{
		return 1;
	}

	// t is supposed to be between 0 and period
	virtual double get_light(size_t __attribute__ ((unused)) day,
		__attribute__ ((unused)) double t) const
	{
		return _level;
	}

	virtual std::string get_name() const
	{
		std::ostringstream oss;
		oss << "ll" << _period << "_" << _level;
		return oss.str();
	}

private:
	double _level;
};


class lightlevelgetter_ldll: public lightlevelgetter
{
public:
	lightlevelgetter_ldll(double period, size_t days_in_ld, double photoperiod,
		size_t days_in_ll, double level):
		lightlevelgetter(period),
		_lightson(.5 * (period - photoperiod)),
		_lightsoff(.5 * (period + photoperiod)),
		_days_in_ld(days_in_ld),
		_days_in_ll(days_in_ll),
		_level(level)
	{
	}

	virtual size_t get_days() const
	{
		return _days_in_ld + _days_in_ll;
	}

	// t is supposed to be between 0 and period
	virtual double get_light(size_t day, double t) const
	{
		if(day < _days_in_ld)
		{
			if(t >= _lightson && t <= _lightsoff)
				return 1000;  // because harvard forest data is about this much.
			return 0;
		}
		else if(day == _days_in_ld && t < _lightson)
			return 0;			// Let the last night of LD run until morning
		return _level;
	}

	virtual std::string get_name() const
	{
		std::ostringstream oss;
		double pp = _lightsoff - _lightson;
		oss << "ld" << pp << "_" << _period - pp << "_ll" << _level;
		return oss.str();
	}

private:
	//light can shift twice each period
	double _lightson;		//light on
	double _lightsoff;	//light off
	size_t _days_in_ld, _days_in_ll;
	double _level;			// level in LL
};



template<class T>
static void use_every_n_element(std::vector<T> &vec, size_t n)
{
	std::vector<T> tmp;
	for(size_t i = 0; i < vec.size(); ++i)
	{
		if(i % n == 0)
			tmp.push_back(vec[i]);
	}
	vec = tmp;
}


class lightlevelgetter_hf: public lightlevelgetter
{
private:
	// Hermite interpolation of function f with derivative fp for x in [0, 1],
	// with function values at 0 and 1 given by f0 and f1, and ditto for fp.
	static inline double hermite_interpolation(double x, double f0, double fp0,
		double f1, double fp1)
	{
		const double xx = x * x;
		return f0 + (f1 - f0) * (3 - 2 * x) * xx +
			(xx - x) * ((fp0 + fp1) * x - fp0);
	}

	// The derivative of the above Hermite interpolation
	static inline double hermite_interpolation_deriv(double x, double f0,
		double fp0, double f1, double fp1)
	{
		const double tmp = 6 * (f1 - f0);
		return fp0 + ((tmp - 4 * fp0 - 2 * fp1) +
			(3 * (fp0 + fp1) - tmp) * x) * x;
	}

	// The second derivative, note: not sure it is scaled correctly, but a
	// quick glance at the results looks right.
	static inline double hermite_interpolation_deriv2(double x, double f0,
		double fp0, double f1, double fp1)
	{
		const double tmp = 6 * (f1 - f0);
		return (tmp - 4 * fp0 - 2 * fp1) + 2 * (3 * (fp0 + fp1) - tmp) * x;
	}

	template<class iterT>
	static double spline_slope(const iterT &i, const iterT &i2)
	{
		return (i2->v - i->v) / (i2->time - i->time);
	}

	// Compute derivative in each point, for a polynomial that will not
	// have any value higher than the highest data point
	template<class vecT>
	void prepare_nonovershooting_spline(vecT &vec)
	{
		assert(vec.size() > 2);
		vec.front().dvdt = spline_slope(vec.begin(), vec.begin() + 1);
		vec.back().dvdt = spline_slope(vec.end() - 1, vec.end() - 2);
		for(typename vecT::iterator i = vec.begin() + 1; i < vec.end() - 1; ++i)
		{
			double k1 = spline_slope(i - 1, i);
			double k2 = spline_slope(i, i + 1);
			if((k1 <= 0 && k2 >=0) || (k1 >= 0 && k2 <=0))
				i->dvdt = 0;
			else
				i->dvdt = 2 / (1 / k1 + 1 / k2);
		}
	}

	struct timecurvepoint_t
	{
		timecurvepoint_t(double t, double _v, double _dvdt):
			time(t), v(_v), dvdt(_dvdt)
		{
		}

		timecurvepoint_t(double t, double _v): time(t), v(_v), dvdt(0)
		{
		}

		inline bool operator<(double t) const { return time < t; }

		double interpolate(const timecurvepoint_t &next, double t) const
		{
			double deltat = next.time - time;
			double tnorm = (t - time) / deltat;
			if(tnorm <= 0)
				return v;
			else if(!(tnorm < 1.))	// Also catches overflows
				return next.v;
			return hermite_interpolation(tnorm, v, dvdt * deltat,
												  next.v, next.dvdt * deltat);
		}

		double interpolate_deriv(const timecurvepoint_t &next, double t) const
		{
			double deltat = next.time - time;
			double tnorm = (t - time) / deltat;
			if(tnorm <= 0)
				return dvdt;
			else if(!(tnorm < 1.))	// Also catches overflows
				return next.dvdt;
			return hermite_interpolation_deriv(tnorm, v, dvdt * deltat,
														  next.v, next.dvdt * deltat) / deltat;
		}

		double interpolate_deriv2(const timecurvepoint_t &next, double t) const
		{
			double deltat = next.time - time;
			double tnorm = (t - time) / deltat;
			if(tnorm <= 0)
				return dvdt;
			else if(!(tnorm < 1.))	// Also catches overflows
				return next.dvdt;
			return hermite_interpolation_deriv2(tnorm, v, dvdt * deltat,
														  next.v, next.dvdt * deltat) / deltat;
		}

		double time, v, dvdt;
	};

	typedef std::vector<timecurvepoint_t> timecurvepointvec_t;

	struct oneday
	{
		oneday() {}

		void push_offset(const timecurvepoint_t& p, double offset)
		{
			data.push_back(timecurvepoint_t(p.time - offset, p.v, p.dvdt));
		}

		timecurvepointvec_t data;
	};

	void read(std::istream &in)
	{
		timecurvepointvec_t data;
		std::string s;
		while(std::getline(in, s))
		{
			if(s.empty() || s[0] == '#')
				continue;
			std::stringstream ss(s);
			double day, sup, sdn;
			int dayi, houri;
			std::string rn;
			ss >> day >> dayi >> houri >> rn >> sup >> sdn;
			data.push_back(timecurvepoint_t((day - 1) * 24, std::max(-sdn, 0.) * _modify_level));
			if(!ss)
				throw std::string("raw data input format error " + s);
		}
		prepare_nonovershooting_spline(data);

		_days.clear();
		double firstt = std::ceil(data.front().time / 24) * 24;
		double nextt = firstt, prevt = -1;
		for(timecurvepointvec_t::const_iterator i = data.begin();
			i != data.end(); ++i)
		{
			double t = i->time;

			// Don't use strange (negative, etc.) points
			if(t < firstt)
				continue;

			// If it is the first point in a new day, then just add it
			if(t < nextt)
			{
				_days.back().push_offset(*i, prevt);
				continue;
			}
			double v = (i - 1)->interpolate(*i, nextt);
			double dvdt = (i - 1)->interpolate_deriv(*i, nextt);
			if(!_days.empty())
				_days.back().push_offset(timecurvepoint_t(nextt, v, dvdt), prevt);
			prevt = nextt;
			nextt += 24;
			_days.push_back(oneday());
			_days.back().push_offset(timecurvepoint_t(prevt, v, dvdt), prevt);
			_days.back().push_offset(*i, prevt);
		}

		// If the last day is not complete, remove it
		if(_days.back().data.back().time != prevt)
			_days.pop_back();
		if(_days.empty())
			throw std::string("realppdata error: no days read");
	}


public:
	lightlevelgetter_hf(const std::string &fname, const std::string &append_str, double mod_intensity):
		lightlevelgetter(24), _append_name(append_str), _modify_level(mod_intensity)
	{
		std::ifstream in(fname);
		if(!in)
			throw "unable to open " + fname;
		read(in);
		use_every_n_element(_days, 2);  // Don't use all data, speed gain
	}

	virtual size_t get_days() const
	{
		return _days.size();
	}


	// t is supposed to be between 0 and period
	double get_light(size_t day, double t) const
	{
		assert(t <= _period && t >= 0);
		timecurvepointvec_t curve = _days.at(day).data;

		timecurvepointvec_t::const_iterator it2 =
			std::lower_bound(curve.begin(), curve.end(), t);
		if(it2 == curve.end())
		{
//			throw bugger(std::string("insanity in ") + __FUNCTION__);
			return curve.back().time;
		}
		if(it2 == curve.begin())
			return it2->v;
		return (it2 - 1)->interpolate(*it2, t);
	}

	virtual std::string get_name() const
	{
		return _append_name;
	}

  private:
	// store time series separated into days
	std::vector<oneday> _days;

	std::string _append_name;

	// factor to multiply light output with
	double _modify_level;
};


#endif
