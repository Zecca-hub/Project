#ifndef LAGEVIEN_SIM
#define LAGEVIEN_SIM

#include <cstdint>
#include <cmath>
#include <vector>

void Force(std::vector<double> *Fx,
           std::vector<double> *Fy,
           std::vector<double> *Fz,
           const std::vector<double> *x,
           const std::vector<double> *y,
           const std::vector<double> *z,
           const double N, const double Z);

double U100(const double x,
            const double y,
            const double z,
            const double Z);

void dU100(double &lap_U100, double &cor_U100,
           const double x, 
           const double y,
           const double z);

double U200(const double x,
            const double y,
            const double z,
            const double Z);

void dU200(double &lap_U200, double &cor_U200,
           const double x, const double y,
           const double z);

void Jastrow(double &f_j, double &f_j_d, double &f_j_dx,
             double &f_j_dy, double &f_j_dz, double &f_j_l,
             const double x_ij, const double y_ij,
             const double z_ij, const double r_ij,
             const double beta);

void PSI2_tot(double &psi, double &psi_b, double &psi_1,
              const std::vector<double> *x,
              const std::vector<double> *y,
              const std::vector<double> *z,
              const std::vector<double> *x1,
              const std::vector<double> *y1,
              const std::vector<double> *z1,
              const double Z, const double N,
              const double Lb, const double Beta,
              const double Beta1);

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
            const double dt, const double N,
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
               const double N, const double dt,
               const double A_r);

void LocalEnergy(double &LE,
                 const std::vector<double> *x,
                 const std::vector<double> *y,
                 const std::vector<double> *z,
                 const double N, const double beta);

#endif
