#include "vector_help.h"

void map(std::vector<double> &vec, const std::function<double(uint32_t)> f)
{
	const uint32_t M{static_cast<uint32_t>(vec.size())};
	for (uint32_t i{0}; i < M; ++i)
	{
		vec[i] = f(i);
	}
}
void arange(std::vector<double> &vec, const double start, const double step)
{
	std::function<double(uint32_t)> f{
		[&](uint32_t i) -> double { return start + i * step; }};
	map(vec, f);
}
void linspace(std::vector<double> &vec, const double xmin, const double xmax)
{
	const double h{(xmax - xmin) / (vec.size() - 1)};
	arange(vec, xmin, h);
}
void fill(std::vector<double> &vec, const double val)
{
	std::function<double(uint32_t)> f{
		[&](uint32_t i) -> double { return val; }};
	map(vec, f);
}

uint32_t ind_min(const std::vector<double> &vec,
				 std::function<double(double)> map)
{

	uint32_t ind{0};
	double val{map ? map(vec[0]) : vec[0]};

	for (uint32_t i{0}; i < static_cast<uint32_t>(vec.size()); ++i)
	{
		if (map != nullptr)
		{
			if (map(vec[i]) < val)
			{
				val = map(vec[i]);
				ind = i;
			}
		}
		else
		{
			if (vec[i] < val)
			{
				val = vec[i];
				ind = i;
			}
		}
	}

	return ind;
}
uint32_t ind_max(const std::vector<double> &vec,
				 std::function<double(double)> map)
{
	if (map != nullptr)
	{
		return ind_min(vec,
					   [&](double x) -> double { return -map(x); });
	}
	else
	{
		return ind_min(vec,
					   [&](double x) -> double { return -x; });
	}
}
