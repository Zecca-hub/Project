#include "lu_dec.h"

LU_Solver::LU_Solver(
	const uint32_t Nx,
	const std::function<double(uint32_t)> qmin,
	const std::function<double(uint32_t)> qmax ):
	Nx{Nx}, qmin{qmin}, qmax{qmax},
	en(Nx-3), e0(Nx-2), ep(Nx-3), xx(Nx-2)
	{
}
LU_Solver::~LU_Solver(){

}

void LU_Solver::step_setup(
	const std::function<double(uint32_t)> dn,
	const std::function<double(uint32_t)> d0,
	const std::function<double(uint32_t)> dp)
	{
	// Calculate the coefficients of the decomposition
	ep[0]=dp(1);
	e0[0]=d0(1);
	en[0]=dn(2)/e0[0];
	for(uint32_t n=1;n<Nx-3;++n){
		ep[n]=dp(n+1);
		e0[n]=d0(n+1)-dn(n+1)*ep[n-1]/e0[n-1];
		en[n]=dn(n+2)/e0[n];
	}
	e0[Nx-3]=d0(Nx-3+1)-dn(Nx-3+1)*ep[Nx-3-1]/e0[Nx-3-1];

	dn_1=dn(1);
	dp_Nx_2=dp(Nx-2);
}

void LU_Solver::step(
	const std::vector<double> *u_p0,
	std::vector<double> *u_p1,
	const uint32_t m0)
	{

	// Eq. 61
	xx[0]=(*u_p0)[1]-dn_1*qmin(m0+1);
	for(uint32_t n=1;n<Nx-3;++n){
		xx[n]=(*u_p0)[n+1]-en[n-1]*xx[n-1];
	}
	xx[Nx-3]=(*u_p0)[Nx-3+1]-dp_Nx_2*qmax(m0+1)-en[Nx-3-1]*xx[Nx-3-1];

	// eq. 64
	(*u_p1)[Nx-1]=qmax(m0+1);
	(*u_p1)[Nx-2]=xx[Nx-3]/e0[Nx-3];
	for(uint32_t n=1;n<Nx-2;++n){
		(*u_p1)[Nx-n-2]=xx[Nx-n-3]/e0[Nx-n-3]
			-ep[Nx-n-3]/e0[Nx-n-3]*(*u_p1)[Nx-n-1];
	}
	(*u_p1)[0]=qmin(m0+1);

}
