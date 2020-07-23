#ifndef __FINDZERO_H__
#define __FINDZERO_H__

#include "differentiate.h"
#include "functions.h"

#include <cstdint>
#include <cmath>
#include <functional>

/* Find zero with Newton Raphson
 * F is the function that this method wants to find the zero of
 * x0 is the starting value for the search
 * epsilon precision i want for the found zero
 * h is the mesh size for the derivation
 */
double findzero_newton_raphson_xeps(std::function<double(double)> F,double x0,double epsilon,double h,
                                    const double xmin=-INFINITY,const double xmax=INFINITY);
double findzero_newton_raphson_yeps(std::function<double(double)> F,double x0,double epsilon,double h);

/* Find zero with bisection
 * F function for the zero
 * xmin and xmax are the boundary
 * epsilon is the accuracy of the zero
 * NOTE: must change sign between xmin and xmax otherwise returns the middle
 */
double findzero_bisection_xeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon);
double findzero_bisection_yeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon);

/* Find zero with secants method
 * F functoin for the zero
 * x0 is the starting point
 * x1 is the other secant starting point
 * epsilon is the precision or digits is the number of significant digits needed
 */
double findzero_secants_xeps(std::function<double(double)> F, double x0,double x1,
    const double epsilon=1e-16,const double xmin=-INFINITY,const double xmax=INFINITY);
double findzero_secants_xdigits(std::function<double(double)> F, double x0,double x1,
    const int digits=5,const double xmin=-INFINITY,const double xmax=INFINITY);


#endif // __FINDZERO_H__