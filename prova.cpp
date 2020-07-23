#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>
#include <random>

void prova(double &x,
           const double y,
           const double z)
{
    x += 3 * y + z;
}

int main()
{
    double a{0};
    for(uint32_t i{0};i<4;i++){

        int y{5};
        int z{5};
        prova(a, y, z);
    }
        std::cout<<a<<std::endl;
        return 0;
}