#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>
#include <random>

#include "langevien_sim.h"

// Evaluation of foundamental state energy of beryllium atom, maybe

int main2()
{

    // Physics constant
    uint32_t N_r{4};                          //# of electron
    uint32_t Z_r{N_r};                        // # of proton
    double e_r{1.6 * 1e-19};                  // C
    double epsilon_0_r{8.8541878128 * 1e-12}; // N/m
    double h_bar_r{1.054571817 * 1e-34};      // Js
    double m_e_r{9.1093837015 * 1e-31};       // Kg
    double R_he_r{112 * 1e-12};               //m
    double V_r{std::pow(2 * R_he_r, 3)};      //m^3
    
    // Conversion value
    double Kx{e_r * e_r * m_e_r / h_bar_r / h_bar_r / 4 / M_PI / epsilon_0_r}; // 1/m, so 1/Kx convert in meter
    double Kf{e_r * e_r * Kx * Kx / 4 / M_PI / epsilon_0_r};                   // N

    // Problem constant
    double dL{R_he_r * Kx}; // atomic "radius"
    double L{6 * dL};       // "lenght" of the simulation box
    double dt{0.001};
    uint32_t M{1000};       // point of simulation

    //todo: find the correct form of A, using the Fokker-Planck equation
    // todo: calculate the variance A of the random force
    // ESTIMATE THIS, maybe it is a constant value 
    double Var{1};          // VarF/Kf/Kf


    // Problem variabiles
    // First position
    std::vector<double> xx0(N_r);
    std::vector<double> *x0=&xx0;
    std::vector<double> yy0(N_r);
    std::vector<double> *y0=&yy0;
    std::vector<double> zz0(N_r);
    std::vector<double> *z0=&zz0;
    // Second position
    std::vector<double> xx1(N_r);
    std::vector<double> *x1=&xx1;
    std::vector<double> yy1(N_r);
    std::vector<double> *y1=&yy1;
    std::vector<double> zz1(N_r);
    std::vector<double> *z1=&zz1;

    // Deterministic Force in position x_n
    std::vector<double> F_xx0(N_r);
    std::vector<double> *F_x0=&F_xx0;
    std::vector<double> F_yy0(N_r);
    std::vector<double> *F_y0=&F_yy0;
    std::vector<double> F_zz0(N_r);
    std::vector<double> *F_z0=&F_zz0;

    // Deterministic Force in position x_n+1
    std::vector<double> F_xx1(N_r);
    std::vector<double> *F_x1=&F_xx1;
    std::vector<double> F_yy1(N_r);
    std::vector<double> *F_y1=&F_yy1;
    std::vector<double> F_zz1(N_r);
    std::vector<double> *F_z1=&F_zz1;

    // Swap puntator
    std::vector<double> *tmp;
    
    // Prabability variabile
    double P;
    double a;

    // Variational parameters
    double beta = 7e-2;
    double beta_1 = 7.1e-2;

    // when we do a cycle over beta, we must nullifier these variabiles
    // Local energy
    double Energy{0};
    // Variance of local energy
    double Energy2;
    double DE{0};
    // gradient of locla energy
    double dE{0};
    double ddE{0};
    double dE_n{0};
    double dE2_n{0};
    double dE_d{0};

    // Prabaility desity function
    double Psi2,Psi2_b,Psi2_1;


    // Particle initzialiazion
    // We put the infinitve massive nucleos in the middle of the box
    // First particle
    (*x0)[0] = 2 * dL;
    (*y0)[0] = -2 * dL;
    (*z0)[0] = -2 * dL;

    // Second particle
    (*x0)[1] = 2 * dL;
    (*y0)[1] = 2 * dL;
    (*z0)[1] = 2 * dL;

    // Third particle
    (*x0)[2] = -2 * dL;
    (*y0)[2] = 2 * dL;
    (*z0)[2] = -2 *dL;

    // Fourth particle
    (*x0)[3] = -2 * dL;
    (*y0)[3] = -2 * dL;
    (*z0)[3] = 2 * dL;


    // Evolution of our system with Langevien mechanics, maybe.

    std::cout << "Lenght conversion factor: " << Kx << '\n'
              << "Force conversion factor: " << Kf << std::endl;


    for (uint32_t i{0}; i < M; i++)
    {

        //todo: a cicle that find the "optimal" dt value

        // Evaluation the initial force
        Force(F_x0, F_y0, F_z0, x0, y0, z0, N_r, Z_r);

        // Evolution of the sistem using the Langevien dynamics
        Langevien(x1, y1, z1, x0, y0, z0, F_x0, F_y0, F_z0, N_r, dt,Var);

        // Evaluation the force
        Force(F_x1, F_y1, F_z1, x1, y1, z1, N_r, Z_r);

        // Estimation of probability for beta, x, beta1, x1
        PSI2_tot(Psi2, Psi2_b, Psi2_1, x0, y0, z0, x1, y1, z1, N_r, Z_r, L, beta, beta_1);//beta[i],beta[i+1]
        

        // Estimation of the acceptance ratio
        P = Prob(x0, y0, z0, x1, y1, z1, F_x0, F_y0, F_z0, F_x1, F_y1, F_z1, Psi2_1, Psi2, dt, N_r, Var);

        a = (double)rand()/RAND_MAX;

        if (a > P)
        {
            // Swap the position for the evolution of system
            tmp = x0;
            x0 = x1;
            x1 = x0;

            tmp = y0;
            y0 = y1;
            y1 = y0;

            tmp = z0;
            z0 = z1;
            z1 = z0;
        }
        // Estimation of probability density
        PSI2_tot(Psi2, Psi2_b, Psi2_1, x0, y0, z0, x1, y1, z1, Z_r, N_r, L, beta, beta_1);

        // Estiamtion of total energy   
        LocalEnergy(Energy,x0,y0,z0,N_r,beta);
        Energy2 += Energy * Energy;

        // Estimation of numerator and denominator for the gradient of E
        dE_n += Psi2_b / Psi2 * Energy;
        dE2_n += Psi2_b / Psi2 * Energy*Energy;
        dE_d += Psi2_b / Psi2;

        // todo: correlation function, so tau
    }

    Energy = Energy / M;
    DE = std::sqrt(1.0 / (M - 1) * (1.0 / M * Energy2 - 1.0 / M / M * Energy * Energy));

    dE = dE_n / dE_d;
    ddE = std::sqrt(1.0 / (M - 1) * (1.0 / M * dE2_n / dE_d - 1.0 / M / M * dE * dE));
    return 0;
}