#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <gsl/gsl_multifit_nlin.h>

int expfunction(const gsl_vector *x, void *params, gsl_vector *f)
{
	const std::vector<double> &data = *(const std::vector<double> *)params;
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	for(size_t i = 0; i < data.size(); ++i)
		gsl_vector_set(f, i, (a * std::exp(b * i)) - data[i]);
	return GSL_SUCCESS;
}

int expderiv(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double a = gsl_vector_get(x, 0);
	double b = gsl_vector_get(x, 1);
	for(size_t i = 0; i < J->size1; ++i)
	{
		double tmp = std::exp(b * i);
		gsl_matrix_set(J, i, 0, tmp);
		gsl_matrix_set(J, i, 1, a * i * tmp);
	}
	return GSL_SUCCESS;
}

int exp_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	expfunction(x, params, f);
	return expderiv(x, params, J);
}

void detrend2(std::vector<double> &data)
{
	// Base initial guess on linear fit
	double sumx = 0;
	for(size_t i = 0; i < data.size(); ++i)
		sumx += data[i];
	
	double parray[] = { sumx / data.size(), 0 };
	gsl_vector_view params = gsl_vector_view_array(parray, 2);

	gsl_multifit_function_fdf fdf;
	fdf.f = &expfunction;
	fdf.df = &expderiv;
	fdf.fdf = &exp_fdf;
	fdf.n = data.size();
	fdf.p = 2;
	fdf.params = &data;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(
		gsl_multifit_fdfsolver_lmsder, data.size(), 2);
	gsl_multifit_fdfsolver_set(solver, &fdf, &params.vector);

	for(int i = 0; i < 100; ++i)
	{
		if(gsl_multifit_fdfsolver_iterate(solver) ||
			gsl_multifit_test_delta(solver->dx, solver->x, 1e-40, 1e-40))
			break;
	}

	double a = gsl_vector_get(solver->x, 0);
	double b = gsl_vector_get(solver->x, 1);
	for(size_t i = 0; i < data.size(); ++i)
		data[i] -= a * std::exp(b * i);

	gsl_multifit_fdfsolver_free(solver);
}


double running_average(const std::vector<double> &data, size_t window_width)
{
	assert(window_width & 1);
	double totsum = 0;
	for(size_t i = window_width / 2; i < data.size() - 1 - window_width / 2; ++i)
		totsum += std::abs(data[i]);

	double deviation = 0;	// from running average
	for(size_t i = 0; i < data.size() - window_width; ++i)
	{
		const double hw = .5 * (window_width - 1);
		double average = 0, wsum = 0;
		for(size_t j = 0; j < window_width; ++j)
		{
			double t = (j - hw) / hw;
//			double w = .25 * (2. - std::abs(t)) * (2. - std::abs(t));
			double w = 1. - std::abs(t);
//			double w = 1.;
			average += w * data[i + j];
			wsum += w;
		}
		average /= wsum;
		double window_center = data[i + window_width/2];
		deviation += std::abs(average - window_center);
	}

	return deviation / totsum;
}

void test(double period, double power)
{
	int pph = 3, win = 48;
	std::vector<double> data((72 + win) * pph);
	for(size_t i = 0; i < data.size(); ++i)
		data[i] = pow(.5 + .5 * sin(i * 2. * M_PI / (period * pph)), power) * exp(i / -120. / pph);

	detrend2(data);
	auto d = running_average(data, win * pph + 1);
	std::cout <<""<< period << "\t" << power << "\t" << d << "\n";
}

int main()
{
	for(double p = 1; p <= 5; p += 1)
	{
		for(double per = 3; per < 100; per += .25)
			test(per, p);
		std::cout << "\n";
	}
}
