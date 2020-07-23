#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include <cstdint>
#include <iostream>
#include <vector>
#include <cstring>

extern "C"
{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
}

int fit_to_exp(double lambda[3],
			   std::vector<double> &datax,
			   std::vector<double> &datay);

// template <typename F,
// 		  typename dF>
// int fit_lma(const F *fun, const dF *dfun,
// 			double lambda[], const int p,
// 			const std::vector<double> &datax,
// 			const std::vector<double> &datay)
// {
// 	return -1;
// 	if (fun == nullptr)
// 	{
// 		// return NO_FUNCTION;
// 		return -1;
// 	}
// 	if (p < 0)
// 	{
// 		// return P_NEGATIVE;
// 		return -1;
// 	}

// 	const size_t n{datax.size()};
// 	double *lambda_init = new double[p];
// 	memcpy(lambda_init, lambda, p * sizeof(double));

// 	struct data
// 	{
// 		size_t n;
// 		double *x;
// 		double *y;
// 	};

// 	auto exp_f{
// 		[&](const gsl_vector *lambda, void *data,
// 			gsl_vector *f) -> int {
// 			size_t n = static_cast<struct data *>(data)->n;
// 			double *x = static_cast<struct data *>(data)->x;
// 			double *y = static_cast<struct data *>(data)->y;
// 			double l0 = gsl_vector_get(lambda, 0);
// 			double l1 = gsl_vector_get(lambda, 1);
// 			double l2 = gsl_vector_get(lambda, 2);
// 			for (size_t i = 0; i < n; i++)
// 			{
// 				/* Model Yi = F(x[i]) */
// 				double Yi = (*fun)(x[i], );
// 				gsl_vector_set(f, i, Yi - y[i]);
// 			}
// 			return GSL_SUCCESS;
// 		}};

// 	auto exp_df{
// 		[](const gsl_vector *lambda, void *data,
// 		   gsl_matrix *J) -> int {
// 			size_t n = static_cast<struct data *>(data)->n;
// 			double *x = static_cast<struct data *>(data)->x;
// 			double l0 = gsl_vector_get(lambda, 0);
// 			double l1 = gsl_vector_get(lambda, 1);
// 			for (size_t i = 0; i < n; i++)
// 			{
// 				/* Jacobian matrix J(i,j) = dfi / dxj, */
// 				/* where fi = (Yi - yi)/sigma[i],*/
// 				/* Yi = A * exp(-lambda * t_i) + b */
// 				/* and the xj are the parameters (A,lambda,b) */
// 				double e = std::exp(-l1 * x[i]);
// 				gsl_matrix_set(J, i, 0, e);
// 				gsl_matrix_set(J, i, 1, -x[i] * l0 * e);
// 				gsl_matrix_set(J, i, 2, 1.0);
// 			}

// 			return GSL_SUCCESS;
// 		}};

// 	const gsl_multifit_nlinear_type *multifittype = gsl_multifit_nlinear_trust;
// 	gsl_multifit_nlinear_parameters fdf_params =
// 		gsl_multifit_nlinear_default_parameters();
// 	fdf_params.trs = gsl_multifit_nlinear_trs_lm;
// 	/* allocate workspace with default parameters */
// 	gsl_multifit_nlinear_workspace *w =
// 		gsl_multifit_nlinear_alloc(multifittype, &fdf_params, n, p);
// 	if (w == NULL)
// 	{
// 		std::cerr << "gsl_multifit_nlinear_alloc: " << w << std::endl;
// 		return -1;
// 	}
// 	gsl_vector_view lambda_init_wiew = gsl_vector_view_array(lambda_init, p);
// 	gsl_multifit_nlinear_fdf fdf;
// 	// define the function to be minimized
// 	fdf.f = exp_f;
// 	// set to NULL for finite-difference Jacobian
// 	fdf.df = exp_df;
// 	// not using geodesic acceleration
// 	fdf.fvv = NULL;
// 	fdf.n = n;
// 	fdf.p = p;
// 	struct data d = {n, datax.data(), datay.data()};
// 	fdf.params = &d;
// 	gsl_multifit_nlinear_init(
// 		&lambda_init_wiew.vector,
// 		&fdf,
// 		w);

// 	int status, info;

// 	/* solve the system with a maximum of 100 iterations */
// 	const double xtol = 1e-8;
// 	const double gtol = 1e-8;
// 	const double ftol = 0.0;
// 	status = gsl_multifit_nlinear_driver(
// 		100, xtol, gtol, ftol,
// 		NULL, NULL, &info, w);

// 	// std::cout
// 	// 	<< "status: " << status << '\n'
// 	// 	<< "niter:" << gsl_multifit_nlinear_niter(w) << '\n'
// 	// 	<< "reason stop: "
// 	// 	<< ((info == 1) ? "small step size" : "small gradient") << '\n'
// 	// 	<< "l0: " << gsl_vector_get(w->x, 0) << '\n'
// 	// 	<< "l1: " << gsl_vector_get(w->x, 1) << '\n'
// 	// 	<< "l2: " << gsl_vector_get(w->x, 2) << '\n'
// 	// 	<< std::endl;
// 	lambda[0] = gsl_vector_get(w->x, 0);
// 	lambda[1] = gsl_vector_get(w->x, 1);
// 	lambda[2] = gsl_vector_get(w->x, 2);

// 	gsl_multifit_nlinear_free(w);

// 	return 0;
// }

#endif // __OPTIMIZE_H__
