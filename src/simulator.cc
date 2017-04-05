#include "simulator.h"

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>


std::ostream *simulator::log = &std::cerr;
void (*simulator::print)(void *params, double t) = 0;

simulator::simulator(size_t dim):
	_dim(dim), _gslstepper(0), _gslctrl(0), _gslevolve(0)
{
	//Runge-Kutta Cash-Karp:
	_gslstepper = gsl_odeiv_step_alloc(gsl_odeiv_step_rkck, dim); 
	
	//Runge-Kutta-Fehlberg:
	//	_gslstepper = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, dim);

	//Bulirsch-Stoer method of Bader and Deuflhard; requires the Jacobian:
	//	_gslstepper = gsl_odeiv_step_alloc(gsl_odeiv_step_bsimp, dim);

	//Keep local error for each step within: eps_abs, eps_rel 
	//(assign temporary values)
	_gslctrl = gsl_odeiv_control_y_new(1e-1, 1e-10);

	//Evolve system a time interval, and keep error according to _gslctrl
	_gslevolve = gsl_odeiv_evolve_alloc(dim);
	assert(_gslstepper && _gslctrl && _gslevolve);
}

//Is simulator is in a vector, which moves to new memory address, is
//it makes a copy. Therefore, transfer state to new copy, and zero old. 
simulator::simulator(const simulator &s):
	_dim(s._dim), _gslstepper(s._gslstepper), _gslctrl(s._gslctrl),
	_gslevolve(s._gslevolve)
{
	s._gslevolve = 0;
	s._gslctrl = 0;
	s._gslstepper = 0;
}

simulator::~simulator()
{
	if(_gslevolve)
		gsl_odeiv_evolve_free(_gslevolve);
	if(_gslctrl)
		gsl_odeiv_control_free(_gslctrl);
	if(_gslstepper)
		gsl_odeiv_step_free(_gslstepper);
}

//steps = iterations, as reference so we can check how many steps were taken
//stepsize = (obvious) adaptive, changes by
//   gsl_*_evolve_apply(_gslevolve, _gslctrl, _gslstepper ...)
//time = starting time
bool simulator::simulate(const control &ctl, const gsl_odeiv_system &sys,
	size_t &steps, double &stepsize, double &time,
	std::vector<double> &integrals) const
{
	assert(sys.dimension == _dim);
	assert(integrals.size() <= _dim);

	gsl_odeiv_control_init(_gslctrl, ctl.eps_abs, ctl.eps_rel, 1, 0);
	gsl_odeiv_step_reset(_gslstepper);
	gsl_odeiv_evolve_reset(_gslevolve);

	if(print)
		print(sys.params, 0);

	std::vector<double> vars(ctl.dataptr, ctl.dataptr + integrals.size());

	double prevprevt = time, prevt = time;
	for(; steps < ctl.maxsteps && time < ctl.maxt; steps++)
	{
		double ss = stepsize;
		//perform the integration: 
		int ret = gsl_odeiv_evolve_apply(_gslevolve, _gslctrl,
			_gslstepper, &sys, &time, ctl.maxt, &ss, ctl.dataptr);
		if(ret != GSL_SUCCESS)
			throw std::string("odeiv failed: ") + gsl_strerror(ret);
		if(time != ctl.maxt || ss > 1e-9)
			stepsize = ss;	// Avoid lowering stepsize very near maxt

		double tw = .5 * (time - prevprevt);
		for(size_t i = 0; i < integrals.size(); ++i)
		{
			integrals[i] += tw * vars[i];
			vars[i] = ctl.dataptr[i];
		}
		prevprevt = prevt;
		prevt = time;

		if(print)
			print(sys.params, time);
	}

	if(time < ctl.maxt)
	{
		if(log)
			*log << "time step limit (" << ctl.maxsteps << ") exceeded at time "
				<<	time << " (target " << ctl.maxt << ")" << std::endl;
		return false;
	}

	double tw = .5 * (time - prevprevt);
	for(size_t i = 0; i < integrals.size(); ++i)
		integrals[i] += tw * vars[i];

	return true;
}


