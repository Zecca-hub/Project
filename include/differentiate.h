#ifndef __DIFFERENTIATE_H__
#define __DIFFERENTIATE_H__

#include <cstdint>
#include <functional>

/* Derivator with 5 points rule
 * of function F
 * around point x0
 * mesh size of h
 */
double derive_5points(std::function<double(double)> F,const double x0,const double h);

#endif // __DIFFERENTIATE_H__