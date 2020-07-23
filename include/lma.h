#include <cstdint>
#include <vector>
#include <cstring>
#include <cmath>

template <typename Vec_x,
		  typename Vec_y,
		  typename F,
		  typename _float>
_float lma_squared_sum(const Vec_x &x,
					   const Vec_y &y,
					   const F &f,
					   const _float *beta)
{
	const size_t N{y.size()};
	_float sum{0}, tmp;
	for (size_t i{0}; i < N; ++i)
	{
		tmp = y[i] - f(x[i], beta);
		sum += tmp * tmp;
	}
	return sum;
}

template <typename F,
		  typename _float>
int lma_gradient_parameters(
	const F &f,
	const _float x0,
	const _float *beta,
	const uint8_t N_params,
	const _float *hs,
	const _float *J)
{
	_float h;
	_float *beta_tmp = new _float[N_params];
	std::memcpy(beta_tmp, beta, N_params * sizeof(_float));

	for (uint8_t i{0}; i < N_params; ++i)
	{
		J[i] = 0;
		h = hs[i];

		// 5 points rule
		beta_tmp[i] = beta[i] - 2 * h;
		J[i] += +f(x0, *beta_tmp);

		beta_tmp[i] = beta[i] - h;
		J[i] += -8 * f(x0, *beta_tmp);

		beta_tmp[i] = beta[i] + h;
		J[i] += +8 * f(x0, *beta_tmp);

		beta_tmp[i] = beta[i] + 2 * h;
		J[i] += -f(x0, *beta_tmp);

		J[i] /= 12 * h;
		beta_tmp[i] = beta[i];
	}
	delete[] beta_tmp;
	return 0;
}

template <typename Vec_x,
		  typename Vec_y,
		  typename F,
		  typename _float>
int lma_curve_fit(
	const Vec_x &data_x,
	const Vec_y &data_y,
	const F &f,
	_float *beta,
	const uint8_t N_params)
{
	const size_t M_data{data_x.size()};

	_float *delta = new _float[N_params];

	_float *J = new _float[M_data * N_params];

	// std::memcpy(tmp, guess, N_params * sizeof(_float));

	_float *hs = new _float[N_params];
	for (uint8_t i{0}; i < N_params; ++i)
	{
		hs[i] = (beta[i] == 0) ? 1e-10 : beta[i] * 1e-10;
	}

	while (true)
	{
		// Construct Jacobian
		for (size_t i{0}; i < M_data; ++i)
		{
			lma_gradient_parameters(
				f, data_x[i], beta, N_params, hs, J + i * N_params);
		}
	}

	delete[] delta, J, tmp, hs;
	return 0;
}
