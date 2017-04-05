#ifndef CLOCKSIM_SIMULATOR_H
#define CLOCKSIM_SIMULATOR_H

#include <iostream>
#include <vector>
#include <cassert>
#include <gsl/gsl_odeiv.h>

class simulator
{
  public:
	struct control
	{
		double eps_abs, eps_rel;
		double maxt;
		size_t maxsteps;
		double *dataptr;
	};

	static std::ostream *log;
	static void (*print)(void *params, double t);

	simulator(size_t dim);
	simulator(const simulator &s);
	virtual ~simulator();
	inline size_t dim() { return _dim; }
	
	// Calls the GSL ODE solver with adaptive step length until time
	// equals ctl.maxt.
	// The number of steps, stepsize and time are updated.
	// The (approximate) integral of the first integrals.size() variables
	// will be added to the corresponding elements in integrals.
	// Returns whether the ODE integration completed successfully
	bool simulate(const control &ctl, const gsl_odeiv_system &sys,
		size_t &steps, double &stepsize, double &time,
		std::vector<double> &integrals) const;

  private:
	size_t _dim;
	mutable gsl_odeiv_step *_gslstepper;
	mutable gsl_odeiv_control *_gslctrl;
	mutable gsl_odeiv_evolve *_gslevolve;
};


#endif
