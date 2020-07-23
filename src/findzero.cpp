#include "findzero.h"

double findzero_newton_raphson_xeps(std::function<double(double)> F,const double x0,const double epsilon,const double h,
									const double xmin,const double xmax){
	double x{x0},corr{F(x)/derive_5points(F,x,h)};

	while( std::abs(corr) > epsilon ){

		corr=F(x)/derive_5points(F,x,h);
		x=x-corr;
		x=(x>=xmax)?xmax:x;
		x=(x<=xmin)?xmin:x;
	}
	return x;
}

double findzero_newton_raphson_yeps(std::function<double(double)> F,const double x0,const double epsilon,const double h){
	double x{x0};

	while( std::abs(F(x)) > epsilon ){

		x=x-F(x)/derive_5points(F,x,h);
	}
	return x;
}

double findzero_bisection_yeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon){
	// must change sign between min and max
	double xmiddle{(xmax+xmin)/2};
	if(F(xmin)*F(xmax) > 0)
		return xmiddle;

	while( std::abs(F(xmiddle)) > epsilon ){
		if (F(xmin)*F(xmiddle) < 0){
			// zero here
			xmax=xmiddle;
		}
		else{
			xmin=xmiddle;
		}
		xmiddle=(xmax+xmin)/2;
	}
	return xmiddle;
}

double findzero_bisection_xeps(std::function<double(double)> F,double xmin,double xmax,const double epsilon){
	// must change sign between min and max
	double xmiddle{(xmax+xmin)/2};
	if(F(xmin)*F(xmax) > 0)
		return xmiddle;

	while( std::abs(xmax-xmin) > epsilon ){
		if (F(xmin)*F(xmiddle) < 0){
			// zero here
			xmax=xmiddle;
		}
		else{
			xmin=xmiddle;
		}
		xmiddle=(xmax+xmin)/2;
	}
	return xmiddle;
}

double findzero_secants_xeps(
	std::function<double(double)> F, double x0,double x1,const double epsilon,
	const double xmin,const double xmax){

	double corr{-F(x1)/( (F(x1)-F(x0))/(x1-x0) )}; 

	while(x1!=x0 && std::abs(corr) > epsilon){
		corr=-F(x1)/( (F(x1)-F(x0))/(x1-x0) );
		x0=x1; 
		x1=x1+corr; 
		x1=x1>xmax?xmax:x1<xmin?xmin:x1;
	}
	return x1;
}

double findzero_secants_xdigits(
	std::function<double(double)> F, double x0,double x1,const int digits,
	const double xmin,const double xmax){
		
	double corr{-F(x1)/( (F(x1)-F(x0))/(x1-x0) )}; 

	while( x1!=x0 && number_of_significant_digits((x1+x0)/2,corr) < digits ){
		corr=-F(x1)/( (F(x1)-F(x0))/(x1-x0) );
		x0=x1; 
		x1=x1+corr; 
		x1=x1>xmax?xmax:x1<xmin?xmin:x1;
	}
	return (x0+x1)/2;
}