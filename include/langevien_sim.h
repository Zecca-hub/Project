#ifndef LAGEVIEN_SIM
#define LAGEVIEN_SIM

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>
#include <random>

void Force(std::vector<double> *Fx,
           std::vector<double> *Fy,
           std::vector<double> *Fz,
           const std::vector<double> *x,
           const std::vector<double> *y,
           const std::vector<double> *z,
           const uint32_t N);

double U100(const double x,
            const double y,
            const double z,
            const uint32_t Z);

void dU100(double &lap_U100, double &cor_U100,
           const double x, 
           const double y,
           const double z, 
           const uint32_t Z);

double U200(const double x,
            const double y,
            const double z,
            const uint32_t Z);

void dU200(double &lap_U200, double &cor_U200,
           const double x, 
           const double y,
           const double z,
           const uint32_t Z);

void Jastrow(double &f_j_d, double &f_j_dx,
             double &f_j_dy, double &f_j_dz, double &f_j_l,
             const double x_ij, const double y_ij,
             const double z_ij, const double r_ij,
             const double beta);

void PSI2_tot(double &psi, double &psi_b,
              double &psi_1, double &psi_1_b,
              const std::vector<double> *x,
              const std::vector<double> *y,
              const std::vector<double> *z,
              const std::vector<double> *x1,
              const std::vector<double> *y1,
              const std::vector<double> *z1,
              const uint32_t N, const double Lb,
              const double Beta, const double Beta1);

double Prob(const std::vector<double> *x0,
            const std::vector<double> *y0,
            const std::vector<double> *z0,
            const std::vector<double> *x1,
            const std::vector<double> *y1,
            const std::vector<double> *z1,
            const std::vector<double> *Fx0,
            const std::vector<double> *Fy0,
            const std::vector<double> *Fz0,
            const std::vector<double> *Fx1,
            const std::vector<double> *Fy1,
            const std::vector<double> *Fz1,
            const double P_1,const double P_0,
            const double dt, const uint32_t N,
            const double A_r);

void Langevien(std::vector<double> *x1,
               std::vector<double> *y1,
               std::vector<double> *z1,
               const std::vector<double> *x0,
               const std::vector<double> *y0,
               const std::vector<double> *z0,
               const std::vector<double> *Fx0,
               const std::vector<double> *Fy0,
               const std::vector<double> *Fz0,
               const uint32_t N, const double dt,
               const double A_r);

void LocalEnergy(double &LE, double &LE1,
                 const std::vector<double> *x,
                 const std::vector<double> *y,
                 const std::vector<double> *z,
                 const uint32_t , const double beta,
                 const double beta1);

#endif
