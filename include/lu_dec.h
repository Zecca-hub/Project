#ifndef __LU_DEC_H__
#define __LU_DEC_H__

#include <vector>
#include <cstdint>
#include <functional>

/* LU decomposition solver:
 * The solution u' has Nx points (i from 0 to Nx-1)
 * of which u[0] and u[Nx-1] are the border points
 * given by qmin(m), qmax(m) functions.
 * The matrix of the single evolution is of the form:
 *  | d0(1) dp(1)                               |  u'[1]        u[1]       dn(0)qmin(m+1)
 *  | dn(2) d0(2) dp(2)                         |  u'[2]        u[2]            0
 *  |        ..............                     | (  ...   ) = ( ...   )-(      0          )
 *  |       dn(Nx-3) d0(Nx-3) dp(Nx-3)          |  u'[Nx-3]     u[Nx-3]         0
 *  |                dn(Nx-2) d0(Nx-2) dp(Nx-2) |  u'[Nx-2]     u[Nx-2]   dp(Nx-1)qmax(m+1)
 * Where u' is the evolved vector (at time m+1)
 * and u is the vector to evolve from (at time m)
 * Info callabilities:
 *  * qmin,qmax are callable for every index m of time you need
 *  * dn(n=[0,...,Nx-2])
 *  * d0(n=[1,...,Nx-2])
 *  * dp(n=[1,...,Nx-1])
 */
class LU_Solver{
	public:
	LU_Solver(
		const uint32_t Nx,
		const std::function<double(uint32_t)> qmin,
		const std::function<double(uint32_t)> qmax);
	~LU_Solver();
	void step_setup(
		const std::function<double(uint32_t)> dn,
		const std::function<double(uint32_t)> d0,
		const std::function<double(uint32_t)> dp
	);
	void step(
		const std::vector<double> *u_p0,
		std::vector<double> *u_p1,
		const uint32_t m0=0
	);

	private:
	const uint32_t Nx;
	std::function<double(uint32_t)> qmin;
	std::function<double(uint32_t)> qmax;
	std::vector<double> en;
	std::vector<double> e0;
	std::vector<double> ep;
	std::vector<double> xx;
	double dn_1;
	double dp_Nx_2;
};

#endif // __LU_DEC_H__
