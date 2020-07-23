#include "differentiate.h"

double derive_5points(std::function<double(double)> F,const double x0,const double h){
	// 5 points derivative
	return ( F(x0-2*h)-8*F(x0-h)+8*F(x0+h)-F(x0+2*h) ) /(12*h);
}