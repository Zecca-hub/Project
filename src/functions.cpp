#include "functions.h"

double randn()
{
	// Statically save the two values
	constexpr double two_pi{M_PI * 2};
	constexpr double scale{1.0 / 32767};

	// If generated==true it already contains the next value
	thread_local double next;
	// Track if already generated
	thread_local bool generated;
	generated = !generated;
	if (generated)
	{
		// Generate new
		double q, w;
		// Don't want a zero, bad for the logarithm
		do
		{
			q = (rand() % 32767) * scale;
		} while (q == 0);
		// Scale the value between [0,1)
		w = (rand() % 32767) * scale;
		// Save for next call
		next = std::sqrt(-2.0 * std::log(q)) * std::cos(two_pi * w);
		// Return value
		return std::sqrt(-2.0 * std::log(q)) * std::sin(two_pi * w);
	}
	else
	{
		// Already good from previous call
		return next;
	}
}

double randu()
{
	constexpr double scale{1.0 / 32767};
	return (rand() % 32767) * scale;
}

int number_of_significant_digits(double x, double corr)
{
	return (int)std::log10(x / corr);
}

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>> &data, std::ofstream &file, const int digits)
{
	uint64_t cols = data.size();

	file << "#";
	// Send column names to the stream
	for (uint64_t j = 0; j < cols; ++j)
	{
		file << data.at(j).first;
		if (j != cols - 1)
			file << ","; // No comma at end of line
	}
	file << "\n";

	for (uint64_t i{0}; i < data[0].second.size(); ++i)
	{
		for (uint64_t j{0}; j < cols; ++j)
		{
			file << boost::format((boost::format("%%.%de") % digits).str()) % data[j].second[i];
			if (j != cols - 1)
				file << ","; // No comma at end of line
		}
		file << '\n';
	}

	return;
}

void tocsv(const std::vector<std::pair<std::string, std::vector<double>>> &data, const std::string &fname, const int digits)
{
	std::ofstream file(fname);
	tocsv(data, file, digits);
	file.close();
}
