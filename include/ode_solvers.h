#ifndef __ODE_SOLVERS_H__
#define __ODE_SOLVERS_H__

#include "vector_help.h"

#include <cstdint>
#include <vector>
#include <functional>

/* Plain implementation of the predictor corrector integrator
 *
 */
double predictor_corrector(std::function<double(double,double)> F,double y,double t,double h);

/* Solve ODE
 * of function F
 * for a step of size h
 * t and y are the initial value
 */
double RK_1_step(std::function<double(double,double)> F, double t, double y, double h);

void RK_2_step(std::function<double(double,double,double)> f1, std::function<double(double,double,double)> f2,
	double t, double x1,double x2, double h,double *x1_,double *x2_);

void RK_2(std::function<double(double,double,double)> f1, std::function<double(double,double,double)> f2,
	double t, double x1,double x2, double h,double *x1_,double *x2_);

#endif // __ODE_SOLVERS_H__
