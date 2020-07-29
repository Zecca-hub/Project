#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>
#include <random>

#include "langevien_sim.h"
#include "vector_help.h"
#include "simulation.h"
#include "optimize.h"

#define URAND() ((double)rand() / ((double)RAND_MAX+1))

// Evaluation of foundamental state energy of beryllium atom, maybe

int main()
{
    // Physics constant
    constexpr uint32_t N_r{4};                          // # of electron
    // uint32_t Z_r{N_r};                        // # of proton
    double e_r{1.6 * 1e-19};                  // C
    double epsilon_0_r{8.8541878128 * 1e-12}; // N/m
    double h_bar_r{1.054571817 * 1e-34};      // Js
    double m_e_r{9.1093837015 * 1e-31};       // Kg
    double R_he_r{112 * 1e-12};               //m
    
    // Conversion value
    double Kx{e_r * e_r * m_e_r / h_bar_r / h_bar_r / 4 / M_PI / epsilon_0_r}; // 1/m
    double Kf{e_r * e_r * Kx * Kx / 4 / M_PI / epsilon_0_r};                   // N

    // Problem constant
    double dL{R_he_r * Kx}; // atomic "radius"
    double L{10 * dL};       // "lenght" of the simulation box
    double dt{1e-5};
    double Var{1};          // given by integrate the langevien equation

	// Indexes sections
    // Number of evolution steps
	constexpr uint32_t M = 200000;
	// Total length of the correlation
	constexpr uint32_t i_correlation_length = 10000;
	// How many means to take prior to apply a correcton to dt
	constexpr uint32_t Delta_correction_n = 100;
    // Thermalization indiex of the system
    uint32_t thermalization_n{0};
    bool thermalized{false};

    // Autocorrelation vector
    std::vector<double> tac(i_correlation_length);
	map(tac, [](uint32_t i) -> double {
		return i;
	});
	
    std::vector<double> Energy_all(M);
    std::vector<double> Energy_(M);
    std::vector<double> dE_n(M);
    std::vector<double> dE_N(M);
    std::vector<double> dE_d(M);
    std::vector<double> dE_D(M);

    std::vector<double> Eac(i_correlation_length);
    std::vector<double> Nac(i_correlation_length);

    // Variabile for dt correction
    double A_mean;
    double P_mean;
    double lambda{0.5};
    std::vector<double> As(Delta_correction_n);
    std::vector<double> Ps(Delta_correction_n);

    // Vector variabiles
    // x_n position
    std::vector<double> xx0(N_r);
    std::vector<double> *x0=&xx0;
    std::vector<double> yy0(N_r);
    std::vector<double> *y0=&yy0;
    std::vector<double> zz0(N_r);
    std::vector<double> *z0=&zz0;
    // x_n+1 position
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

    // Varaitional parameters vector and energy vector
    // Dimension
    double d{35};
    // Parameter vector
    std::vector<double> beta(d);
    arange(beta,0.07,0.002);
    // linspace(beta, 0.007, 0.13);

    // double beta1{0.0012};
    // double beta2{0.0014};

    // Prepare energy output
    // File for autocorrelation
    std::ofstream ENERGY_FILE_AC("data/Eac.dat");
    ENERGY_FILE_AC << "step Eac";
    // File for energy 
    std::ofstream ENERGY_FILE("data/Energy.dat");
    ENERGY_FILE << "beta E dE";
    // File for gradient of energy
    std::ofstream GENERGY_FILE("data/Genergy.dat");
    GENERGY_FILE << "beta GE dGE";

    std::cout << '\n'
              << " Parameter of simulation" << '\n'
              << "Lenght conversion factor: " << Kx << '\n'
              << "Force conversion factor: " << Kf << '\n'
              << "Lenght of box: " << L << '\n'
              << "Evolution step lenght: " << dt << std::endl;

    // Evolution of our system with Langevien mechanics
    for (uint32_t j{0}; j < beta.size()-1; j++)
    {
        std::cout << '\n'
                  << '\n'
                  << '\n'
                  << "###" << '\t' << "Simulation for beta value: " << beta[j] << '\t' << " ###" << '\n'
                  << "Progress for beta " << j + 1 << " of " << d << std::endl;

        // Particle initzialiazion
        // We put the infinitve massive nucleos in the middle of the box
        // First particle
        (*x0)[0] = 2.1 * dL;
        (*y0)[0] = -2 * dL;
        (*z0)[0] = -2 * dL;

        // Second particle
        (*x0)[1] = 2  * dL;
        (*y0)[1] = 2.5 * dL;
        (*z0)[1] = 2 * dL;

        // Third particle
        (*x0)[2] = -2.8 * dL;
        (*y0)[2] = 2 * dL;
        (*z0)[2] = -2.6 * dL;

        // Fourth particle
        (*x0)[3] = -2 * dL;
        (*y0)[3] = -2.7 * dL;
        (*z0)[3] = 2 * dL;

        double Energy{0};
        double Energy_b{0};

        // for x0
        double Psi2{1};
        double Psi2_b{1};

        // for x1
        double Psi2_1{1};
        double Psi2_1_b{1};

        for (uint32_t i{1}; i < M; i++)
        {
            // Evaluation the initial force
            Force(F_x0, F_y0, F_z0, x0, y0, z0, N_r);

            // Evolution of the sistem using the Langevien dynamics
            Langevien(x1, y1, z1, x0, y0, z0, F_x0, F_y0, F_z0, N_r, dt, Var);

            // Evaluation of the force
            Force(F_x1, F_y1, F_z1, x1, y1, z1, N_r);

            // Estimation of probability for x0 (beta and beta') and x1 (beta and beta')
            PSI2_tot(Psi2, Psi2_b, Psi2_1, Psi2_1_b, x0, y0, z0, x1, y1, z1, N_r, L, beta[j], beta[j+1]);

            // Estimation of the acceptance ratio
            P = Prob(x0, y0, z0, x1, y1, z1, F_x0, F_y0, F_z0, F_x1, F_y1, F_z1, Psi2_1, Psi2, dt, N_r, Var);

            As[i % Delta_correction_n] = P;
                    

            // Corrections to dt
            if (!thermalized)
            {
                if (i % Delta_correction_n == 0)
                {
                    A_mean = mean(As);
                    if (A_mean > 0.9 && A_mean < 1)
                    {
                        thermalized = true;
                        thermalization_n = i;
                    }
                    else
                    {
                        // Correction
                        dt *= 1 + lambda * (A_mean - 0.5);
                    }
                }
            }

            a = URAND();

            if (a < P)
            {
                // Swap the position for the evolution of system
                tmp = x0;
                x0 = x1;
                x1 = tmp;

                tmp = y0;
                y0 = y1;
                y1 = tmp;

                tmp = z0;
                z0 = z1;
                z1 = tmp;

                // Swap the probability
                Psi2 = Psi2_1;
                Psi2_b = Psi2_1_b;
            }

            Ps[i % Delta_correction_n] = Psi2_1_b / Psi2_1;

            if (i >= thermalization_n)
            {
                // Estiamtion of total energy
                LocalEnergy(Energy, Energy_b, x0, y0, z0, N_r, beta[j], beta[j+1]);
                Energy_all[i] = Energy;

                // Estimation of numerator and denominator for the gradient of E todo.
                dE_n[i] = Psi2_b / Psi2 * Energy_b;
                dE_d[i] = Psi2_b / Psi2;
            }

            //Print progress
            if ((i % 10000) == 0)
            {
                A_mean = mean(As);
                P_mean = mean(Ps);
                
                std::cout
                    << "Step: " << i
                    << "\tP: " << (i * 100) / M << "%"
                    << "\tE: " << Energy
                    << "\tNum: " << dE_n[i]
                    << "\tDen: " << dE_d[i]
                    << "\tA_mean: " << A_mean
                    << "\tP_mean: " << P_mean
                    << std::endl;
            }
        }

        // Resize of vecotr
        Energy_.resize(M - thermalization_n);
        dE_N.resize(M - thermalization_n);
        dE_D.resize(M - thermalization_n);

        for (uint32_t i{0}; i < Energy_.size(); ++i)
        {
            Energy_[i] = Energy_all[i + thermalization_n];
            dE_N[i] = dE_n[i + thermalization_n];
            dE_D[i] = dE_d[i + thermalization_n];
        }

        // Autocorrelation of the local energy
        autocorrelation(Eac, Energy_);
        autocorrelation(Nac, dE_N);

        for (uint32_t k{0}; k < tac.size(); ++k)
        {
            ENERGY_FILE_AC << tac[k] << " " << Eac[k] << '\n';
        }

        // Estimation of tau with an exponential fit
        double par_E[1] = {(Eac[0] - Eac[1000]) / 1000};
        fit_to_exp_1par(par_E, tac, Eac);

        double par_N[1] = {(Nac[0] - Nac[1000]) / 1000};
        fit_to_exp_1par(par_N, tac, Nac);
        
        // Number of dependent points
        double tau_E = (1.0 / par_E[0]);
        double tau_N = (1.0 / par_N[0]);

        // Mean value with its variance of local energy
        double mean_E = mean(Energy_);
        double Dmean_E = tau_E / Energy_.size() * variance(Energy_, 1);

        std::cout << "Mean value of energy: " << mean_E << "+/-" << std::sqrt(Dmean_E) << std::endl;

        ENERGY_FILE << beta[j] << " " << std::fabs(mean_E) << " " << std::sqrt(Dmean_E) << '\n';

        // // Mean value with its variance of gradient in parameter of local energy
        double mean_E_beta = mean(dE_N) / mean(dE_D);
        double Dmean_E_beta = tau_N / dE_N.size() * variance(dE_N) / mean(dE_D) / mean(dE_D);

        double mean_GE = (mean_E_beta - mean_E) / (beta[j+1] - beta[j]);
        double Dmean_GE = (Dmean_E + Dmean_E_beta) / (beta[j+1] - beta[j]) / (beta[j+1] - beta[j]);

        std::cout << "Mean value of gradient of energy: " << mean_GE << "+/-" << std::sqrt(Dmean_GE) << std::endl;

        GENERGY_FILE << beta[j] << " " << mean_GE << " " << std::sqrt(Dmean_GE) << '\n';
    }

    ENERGY_FILE.close();

    return 0;
}