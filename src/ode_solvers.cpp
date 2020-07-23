#include "ode_solvers.h"

double predictor_corrector(std::function<double(double,double)> F,double y,double t,double h){
	double k1{ F(t,y) };
	double y_pred{ y+h*k1 };
	double k2{ F(t+h,y_pred) };
	return y + h*0.5*(k1+k2);
}

double RK_1_step(const std::function<double(double,double)> F, const double t, const double y, const double h){

	double k1{ F(t,y) };
	double k2{ F(t+h/2,y+k1/2) };
	double k3{ F(t+h/2,y+k2/2) };
	double k4{ F(t+h,y+k3) };

	return y+h/6*(k1+2*k2+2*k3+k4);
}

void RK_2_step(const std::function<double(double,double,double)> f1,
	const std::function<double(double,double,double)> f2,
	const double t, const double x1, const double x2,
	const double h, double *x1_, double *x2_){

	double k11{f1(t,x1,x2)};
	double k21{f2(t,x1,x2)};

	double k12{f1(t+h/2,x1+h/2*k11,x2+h/2*k21)};
	double k22{f2(t+h/2,x1+h/2*k11,x2+h/2*k21)};

	double k13{f1(t+h/2,x1+h/2*k12,x2+h/2*k22)};
	double k23{f2(t+h/2,x1+h/2*k12,x2+h/2*k22)};

	double k14{f1(t+h,x1+h*k13,x2+h*k23)};
	double k24{f2(t+h,x1+h*k13,x2+h*k23)};

	*x1_=x1+h/6*(k11+2*k12+2*k13+k14);
	*x2_=x2+h/6*(k21+2*k22+2*k23+k24);
}

void RK_2(const std::function<double(double,double,double)> f1,
	const std::function<double(double,double,double)> f2,
	const double start_t, const double start_x1, const double start_x2,
	const double stop_t, std::vector<double> &x1, std::vector<double> &x2,
	const uint32_t steps=10){

	x1.clear();
	x1.resize(steps);
	x2.clear();
	x2.resize(steps);

	std::vector<double> t(steps);
	double h{(stop_t-start_t)/steps};

	arange(t,start_t,h);
	x1[0]=start_x1;
	x2[0]=start_x2;

	for (uint32_t i{0};i<steps-1;++i){
		RK_2_step(f1,f2,t[i],x1[i],x2[i],h,&(x1[i+1]),&(x2[i+1]));
	}
}

// void runge_kutta_vector(uint8_t n, const std::function<double(double,const double*)> *f,
// 						double t,const double *x,double h,const double *x_){
// 	double *k1{ new double[n] };
// 	double *k2{ new double[n] };
// 	double *k3{ new double[n] };
// 	double *k4{ new double[n] };

// 	for (uint8_t i{0};i<n;++i){
// 		k1[i]=(f[i])(t,x);
// 		k2[i]=(f[i])(t+h/2,x+h/2*k1[i])
// 	}
// }
