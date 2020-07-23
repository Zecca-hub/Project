#include "langevien_sim.h"

// Evalute the force of the system
void Force(std::vector<double> *Fx,
           std::vector<double> *Fy,
           std::vector<double> *Fz,
           const std::vector<double> *x,
           const std::vector<double> *y,
           const std::vector<double> *z,
           const double N, const double Z)
{
    
    double xx;
    double yy;
    double zz;
    double r;

    for (uint32_t i = 0; i < N; i++)
    {
        (*Fx)[i] = 0;
        (*Fy)[i] = 0;
        (*Fz)[i] = 0;

        r = std::sqrt((*x)[i] * (*x)[i] + (*y)[i] * (*y)[i] + (*z)[i] * (*z)[i]);

        // Coulomb force
        (*Fx)[i] += Z / r / r * (*x)[i] / r;
        (*Fy)[i] += Z / r / r * (*y)[i] / r;
        (*Fz)[i] += Z / r / r * (*z)[i] / r;

        for (uint32_t j = 0; j < N; j++)
        {
            if(j!=i)
            {            
                // Evaluate the relative position
                xx = (*x)[i] - (*x)[j];
                yy = (*y)[i] - (*y)[j];
                zz = (*z)[i] - (*z)[j];

                r = std::sqrt(xx * xx + yy * yy + zz * zz);
                
                // Relative coulomb force(i put the factor 1/2 because the sum is over all i and j)
                (*Fx)[i] += -1.0 / 2 / r / r * xx / r;
                (*Fy)[i] += -1.0 / 2 / r / r * yy / r;
                (*Fz)[i] += -1.0 / 2 / r / r * zz / r;
            }
        }
    }
}

// Different part of the wavefunction
// Wavefunction of isolated electron in 100 state
double U100(const double x,
            const double y,
            const double z,
            const double Z)
{
    double r;
    double psi;

    r = std::sqrt(x * x + y * y + z * z);
    psi = 1.0 / std::sqrt(M_PI) * std::exp(-r * Z);

    return psi;
}

// Derivation of U100 wavefunction
void dU100(double &lap_U100, double &cor_U100,
           const double x, 
           const double y,
           const double z)
{
    // lap_U100 is the laplacian of U100
    // cor_U100 is the sum of the differenzation of U100
    double r;
    
    r = std::sqrt(x * x + y * y + z * z);

    lap_U100 = 8.0 / std::sqrt(M_PI) * std::exp(-4 * r);
    cor_U100 = -4.0 / std::sqrt(M_PI) * (x + y + z) / r * std::exp(-4 * r);
}

// Wavefunction of isolated electron in 200 state
double U200(const double x,
            const double y,
            const double z,
            const double Z)
{
    double r;
    double psi;

    r = std::sqrt(x * x + y * y + z * z);
    psi = 1.0 / 4.0 / std::sqrt(2 * M_PI) * (2 - r) * std::exp(-r * Z / 2);

    return psi;
}

// Derivation of U200 wavefunction
void dU200(double &lap_U200, double &cor_U200,
           const double x, const double y,
           const double z)
{
    // lap_U100 is the laplacian of U100
    // cor_U100 is the sum of the differenzation of U100
    double r;
    
    r = std::sqrt(x * x + y * y + z * z);

    lap_U200 = 2.0 / 3 / std::sqrt(2 * M_PI) * std::exp(-2 * r) * (1 - 1.0 / r);
    cor_U200 = -4.0 / std::sqrt(M_PI) * (x + y + z) / r * std::exp(-4 * r);
}

void Jastrow(double &f_j, double &f_j_d, double &f_j_dx,
             double &f_j_dy, double &f_j_dz, double &f_j_l,
             const double x_ij, const double y_ij,
             const double z_ij, const double r_ij,
             const double beta)
{
    // Jastrow function
    f_j += r_ij / (1 + beta * r_ij);

    // Sum of first derivation of Jastrow function
    f_j_d += (x_ij + y_ij + z_ij) / (r_ij + beta * r_ij) /
            (r_ij + beta * r_ij) / r_ij;

    // Derivation of Jastrow in x,y and z
    f_j_dx += x_ij / (1 + beta * r_ij) / (1 + beta * r_ij) / r_ij;
    f_j_dy += y_ij / (1 + beta * r_ij) / (1 + beta * r_ij) / r_ij;
    f_j_dz += z_ij / (1 + beta * r_ij) / (1 + beta * r_ij) / r_ij;

    // Laplacian of Jastrow function
    f_j_l += (3 * r_ij - r_ij * r_ij * ((1 + beta * r_ij) / r_ij + 2 * beta)) /
            ((1 + beta * r_ij) * (1 + beta * r_ij) * (1 + beta * r_ij) * r_ij * r_ij);
}

// PDF
// change this one, implementation of P(x1)
void PSI2_tot(double &psi, double &psi_b, double &psi_1,
              const std::vector<double> *x,
              const std::vector<double> *y,
              const std::vector<double> *z,
              const std::vector<double> *x1,
              const std::vector<double> *y1,
              const std::vector<double> *z1,
              const double Z, const double N,
              const double Lb, const double Beta,
              const double Beta1)
{
    double xx;
    double x_ij;
    double yy;
    double y_ij;
    double zz;
    double z_ij;

    // For two different beta
    double psi1{1};
    double psi2{1};

    double exp_j{0};
    double exp_j_b{0};

    double A0, B0, C0, D0, E0, F0, G0, H0;

    double r_ij;

    // For the different position
    double psi11{1};
    
    double A1, B1, C1, D1, E1, F1, G1, H1;

    for (uint32_t i{0}; i < N; i++)
    {

        // Control if a particle is go out from the box side Lb
        if ((*x)[i] > Lb || (*y)[i] > Lb || (*z)[i] > Lb)
        {
            psi1 *= 0;
            psi2 *= 0;
            psi11 *=0;
        }

        // Estimetion of Jastrow function
        for (uint32_t j{0}; j < N; j++)
        {
            if (i != j)
            {

                x_ij = (*x)[i] - (*x)[j];
                y_ij = (*y)[i] - (*y)[j];
                z_ij = (*z)[i] - (*z)[j];

                r_ij = std::sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij);

                exp_j += 1.0 / 2 * r_ij / (1 + Beta * r_ij);
                exp_j_b += 1.0 / 2 * r_ij / (1 + Beta1 * r_ij);

                x_ij = (*x1)[i] - (*x1)[j];
                y_ij = (*y1)[i] - (*y1)[j];
                z_ij = (*z1)[i] - (*z1)[j];

                r_ij = std::sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij);

                exp_j += 1.0 / 2 * r_ij / (1 + Beta * r_ij);
            }
        }

        psi1 *= std::exp(exp_j);
        psi2 *= std::exp(exp_j_b);

        psi11 *= std::exp(exp_j);

        if (i == 0)
        {
            // For x0
            A0 = U100((*x)[i], (*y)[i], (*z)[i], Z);
            B0 = U200((*x)[i], (*y)[i], (*z)[i], Z);

            // For x1
            A1 = U100((*x1)[i], (*y1)[i], (*z1)[i], Z);
            B1 = U200((*x1)[i], (*y1)[i], (*z1)[i], Z);
        }

        if (i == 1)
        {
            // For x0
            C0 = U100((*x)[i], (*y)[i], (*z)[i], Z);
            D0 = U200((*x)[i], (*y)[i], (*z)[i], Z);

            // For x1
            C1 = U100((*x1)[i], (*y1)[i], (*z1)[i], Z);
            D1 = U200((*x1)[i], (*y1)[i], (*z1)[i], Z);
        }

        if (i == 2)
        {
            // For x0
            E0 = U100((*x)[i], (*y)[i], (*z)[i], Z);
            F0 = U200((*x)[i], (*y)[i], (*z)[i], Z);

            // For x1
            E1 = U100((*x1)[i], (*y1)[i], (*z1)[i], Z);
            F1 = U200((*x1)[i], (*y1)[i], (*z1)[i], Z);
        }

        if (i == 3)
        {
            // For x0
            G0 = U100((*x)[i], (*y)[i], (*z)[i], Z);
            H0 = U200((*x)[i], (*y)[i], (*z)[i], Z);

            // For x1
            G1 = U100((*x1)[i], (*y1)[i], (*z1)[i], Z);
            H1 = U200((*x1)[i], (*y1)[i], (*z1)[i], Z);
        }
    }
    // For x0
    psi1 *= (A0 * D0 - B0 * C0) * (E0 * H0 - F0 * G0);
    psi2 *= (A0 * D0 - B0 * C0) * (E0 * H0 - F0 * G0);

    // For x1
    psi11 *= (A1 * D1 - B1 * C1) * (E1 * H1 - F1 * G1);

    // Calcolous of probability density punction for 2 different value of Beta
    psi = psi1 * psi1;
    psi_b = psi2 * psi2;

    // Calcolous of probability density function for x1
    psi_1 = psi11*psi11;
}

// Acceptance rate
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
            const double A_r)
{
    // i0 is the old position (with i =x,y,z)
    // i1 is the new position (with i =x,y,z)
    // Fi0 is the old force (with i =x,y,z)
    // Fi1 is the new force (with i =x,y,z)
    // P0 is the probability in x_n
    // P1 is the probability in x_n+1

    double Nu;
    
    double P;
    double A;
    double T;
    
    double exp_x;
    double exp_y;
    double exp_z;
   
        for (uint32_t i{0}; i < N; i++)
    {
        // Evaluation of exponent of transition matrix
        exp_x = dt * dt * ((*Fx0)[i] * (*Fx0)[i] - (*Fx1)[i] * (*Fx1)[i]) 
                + 2 * dt * ((*x0)[i] - (*x1)[i]) * ((*Fx0)[i] + (*Fx1)[i]);
        exp_y = dt * dt * ((*Fy0)[i] * (*Fy0)[i] - (*Fy1)[i] * (*Fy1)[i]) 
                + 2 * dt * ((*y0)[i] - (*y1)[i]) * ((*Fy0)[i] + (*Fy1)[i]);
        exp_z = dt * dt * ((*Fz0)[i] * (*Fz0)[i] - (*Fz1)[i] * (*Fz1)[i]) 
                + 2 * dt * ((*z0)[i] - (*z1)[i]) * ((*Fz0)[i] + (*Fz1)[i]);
    }

    Nu = (exp_x + exp_y + exp_z) / (4 * dt / A_r);

    // Estiamtion of acceptance rate
    // put the value of P(x1) and P(x0)
    A = P_1 / P_0 * std::exp(Nu);

    // Probability 
    P = std::min(A, 1.0);

    return P;
}


// Evolution of the system with Langevien dynamics
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
               const double A_r)
{

    double a, b, c, d, A, B, C;
    for (uint32_t j{0}; j < N; j++)
    {
        // Box-Muller equation for the random force
        a = (double)rand() / RAND_MAX;
        b = (double)rand() / RAND_MAX;
        c = (double)rand() / RAND_MAX;
        d = (double)rand() / RAND_MAX;

        // A,B,C have a normal distrubution with variance 2*A*dt
        A = std::sqrt(-4 * A_r * dt * std::log10(1 - a)) * std::cos(2 * M_PI * b);
        B = std::sqrt(-4 * A_r * dt * std::log10(1 - a)) * std::sin(2 * M_PI * b);
        C = std::sqrt(-4 * A_r * dt * std::log10(1 - c)) * std::cos(2 * M_PI * d);

        // Evolution of the sistem
        (*x1)[j] = (*x0)[j] + dt * (*Fx0)[j] + dt * A;
        (*y1)[j] = (*y0)[j] + dt * (*Fy0)[j] + dt * B;
        (*z1)[j] = (*z0)[j] + dt * (*Fz0)[j] + dt * C;
    }
}

// Evaluate the local energy of the sistem
void LocalEnergy(double &LE,
                 const std::vector<double> *x,
                 const std::vector<double> *y,
                 const std::vector<double> *z,
                 const double N, const double beta)
{
    double LE{0};
    double r;
    double xx, yy, zz;

    double f_j, f_j_d, f_j_dx, f_j_dy, f_j_dz, f_j_l;

    double U_lap10, U_cor10, U_lap20, U_cor20;
    double U_lap11, U_cor11, U_lap21, U_cor21;
    double U_lap12, U_cor12, U_lap22, U_cor22;
    double U_lap13, U_cor13, U_lap23, U_cor23;

    double T0, T1, T2, T3;
    double D0, D1;

    double A, B, C, D, E, F, G, H;
    double a, b, c, d, e, f, g, h, i, j, l, m, n, o, p, q, s, t, u, v;

    for (uint32_t i{0}; i < N; i++)
    {
        //Potential term of local energy
        r = std::sqrt((*x)[i] * (*x)[i] + (*y)[i] * (*y)[i] + (*z)[i] * (*z)[i]);

        LE += N / r;

        for (uint32_t j{0}; j < N; j++)
        {
            if (j != i)
            {
                //Interacation term electron-electron
                xx = (*x)[i] - (*x)[j];
                yy = (*y)[i] - (*y)[j];
                zz = (*z)[i] - (*z)[j];

                r = std::sqrt(xx * xx + yy * yy + zz * zz);

                // factor 0.5 for the sum
                LE += N / 2 / r;

                // Calcolus of Jastrow factor and its derivate
                Jastrow(f_j, f_j_d, f_j_dx, f_j_dy, f_j_dz, f_j_l, xx, yy, zz, r, beta);
                // todo: this don't sum the term
            }
        }

        // Term to estimate che kinetic part of local energy
        if (i == 0)
        {
            dU100(U_lap10, U_cor10, (*x)[i], (*y)[i], (*z)[i]);
            dU200(U_lap20, U_cor20, (*x)[i], (*y)[i], (*z)[i]);
            A = U100((*x)[i], (*y)[i], (*z)[i], 4);
            B = U200((*x)[i], (*y)[i], (*z)[i], 4);
            a = f_j_d;
            b = f_j_dx;
            c = f_j_dy;
            d = f_j_dz;
            e = f_j_l;
        }
        if (i == 1)
        {
            dU100(U_lap11, U_cor11, (*x)[i], (*y)[i], (*z)[i]);
            dU200(U_lap21, U_cor21, (*x)[i], (*y)[i], (*z)[i]);
            C = U100((*x)[i], (*y)[i], (*z)[i], 4);
            D = U200((*x)[i], (*y)[i], (*z)[i], 4);
            f = f_j_d;
            h = f_j_dx;
            g = f_j_dy;
            i = f_j_dz;
            j = f_j_l;
        }
        if (i == 2)
        {
            dU100(U_lap12, U_cor12, (*x)[i], (*y)[i], (*z)[i]);
            dU200(U_lap22, U_cor22, (*x)[i], (*y)[i], (*z)[i]);
            E = U100((*x)[i], (*y)[i], (*z)[i], 4);
            F = U200((*x)[i], (*y)[i], (*z)[i], 4);
            l = f_j_d;
            m = f_j_dx;
            n = f_j_dy;
            o = f_j_dz;
            p = f_j_l;
        }
        if (i == 3)
        {
            dU100(U_lap13, U_cor13, (*x)[i], (*y)[i], (*z)[i]);
            dU200(U_lap23, U_cor23, (*x)[i], (*y)[i], (*z)[i]);
            G = U100((*x)[i], (*y)[i], (*z)[i], 4);
            H = U200((*x)[i], (*y)[i], (*z)[i], 4);
            q = f_j_d;
            v = f_j_dx;
            s = f_j_dy;
            t = f_j_dz;
            u = f_j_l; 
        }
    }

    D0 = A * F - E * B;
    D1 = C * H - D * G;
    
    // Term of kinetic energy of particle 1
    T0 = (U_lap10 * F - U_lap20 * E + 2 * (U_cor10 * F - E * U_cor20) * a) / D0 + e + b * b + c * c + d * d;
    // Term of kinetic energy of particle 2
    T1 = (U_lap11 * H - U_lap21 * G + 2 * (U_cor11 * H - G * U_cor12) * f) / D1 + j + h * h + g * g + i * i;
    // Term of kinetic energy of particle 3
    T2 = (U_lap22 * A - U_lap12 * B + 2 * (U_cor22 * A - U_cor21 * B) * l) / D0 + p + m * m + n * n + o * o;
    // Term of kinetic energy of particle 4
    T3 = (U_lap23 * C - U_lap13 * D + 2 * (U_cor23 * C - U_cor13 * D) * q) / D1 + u + v * v + s * s + t * t;

    LE += T0 + T1 + T2 + T3;
}

