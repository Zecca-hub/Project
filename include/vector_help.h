#ifndef __VECTOR_HELP_H__
#define __VECTOR_HELP_H__

#include <cstdint>
#include <vector>
#include <functional>

// vector manipulation
void map(std::vector<double> &vec, std::function<double(uint32_t)> f);
void arange(std::vector<double> &vec, const double start, const double step);
void linspace(std::vector<double> &vec, const double xmin, const double xmax);
void fill(std::vector<double> &vec, double val);

uint32_t ind_min(const std::vector<double> &vec,
				 std::function<double(double)> map = nullptr);
uint32_t ind_max(const std::vector<double> &vec,
				 std::function<double(double)> map = nullptr);

template <typename T>
T mean(const std::vector<T> &x)
{
	T r{0};
	for (const T &v : x)
	{
		r += v;
	}
	return r / x.size();
}

template <typename T>
T variance(const std::vector<T> &x, const int ddof = 0)
{
	T s1{0};
	T s2{0};
	for (const T &v : x)
	{
		s1 += v;
		s2 += v * v;
	}
	return s2 / (x.size() - ddof) - s1 * s1 / (x.size() * (x.size() - ddof));
}
#endif // __VECTOR_HELP_H__
