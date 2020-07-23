#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#include <cstdint>
#include <vector>
#include <functional>

/* Integrate with simpson cubic rule
 * the function F
 * from xmin to xmax
 * with M steps
 */
double integrator_simpson_cubic(std::function<double(double)> F,const double xmin,const double xmax,const uint32_t M);

#endif // __INTEGRATE_H__