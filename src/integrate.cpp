#include "integrate.h"

double integrator_simpson_cubic(std::function<double(double)> F,const double xmin,const double xmax,const uint32_t M){
	// std::vector<double> f(M+1);
	// std::vector<double> x(M);
	double h{(xmax-xmin)/M};

	// for (uint32_t i{0}; i<=M; ++i){
	// 	f[i]=F(xmin+i*h);
	// }

	double result{0};
	// result+=f[0]+f[M];
	result+=F(xmin)+F(xmax);
	for (uint32_t i{1}; i<=M/2-1; ++i){
		// result+=4*f[2*i-1]+2*f[2*i];
		result+=4*F(xmin+2*h*i-h)+2*F(2*h*i);
	}
	// result+=4*f[M-1];
	result+=4*F(xmax-h);

	return h/3*result;
}