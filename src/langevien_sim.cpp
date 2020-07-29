#include "langevien_sim.h"

// Evalute the force of the system
void Force(std::vector<double> *Fx,
           std::vector<double> *Fy,
           std::vector<double> *Fz,
           const std::vector<double> *x,
           const std::vector<double> *y,
           const std::vector<double> *z,
           const uint32_t N)
{
    uint32_t Z{N};
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
        (*Fx)[i] += -(Z / r / r * (*x)[i] / r);
        (*Fy)[i] += -(Z / r / r * (*y)[i] / r);
        (*Fz)[i] += -(Z / r / r * (*z)[i] / r);

        for (uint32_t j = 0; j < N; j++)
        {
            if (i != j)
            {
                // Evaluate the relative position
                xx = (*x)[i] - (*x)[j];
                yy = (*y)[i] - (*y)[j];
                zz = (*z)[i] - (*z)[j];

                r = std::sqrt(xx * xx + yy * yy + zz * zz);

                // Relative coulomb force
                (*Fx)[i] += 1.0 / r / r * xx / r;
                (*Fy)[i] += 1.0 / r / r * yy / r;
                (*Fz)[i] += 1.0 / r / r * zz / r;
            }
        }
    }
}

// Different part of the wavefunction
// Wavefunction of isolated electron in 100 state
double U100(const double x,
            const double y,
            const double z,
            const uint32_t Z)
{
    double r;
    double psi;

    r = std::sqrt(x * x + y * y + z * z);
    psi = std::pow(Z, 3.0 / 2) / std::sqrt(M_PI) * std::exp(-r * Z);

    return psi;
}

// Derivation of U100 wavefunction
void dU100(double &lap_U100, double &cor_U100,
           const double x, 
           const double y,
           const double z, 
           const uint32_t Z)
{
    // lap_U100 is the laplacian of U100
    // cor_U100 is the sum of the differenzation of U100
    double r;
    
    r = std::sqrt(x * x + y * y + z * z);

    lap_U100 = std::pow(Z, 5.0 / 2) / std::sqrt(M_PI) * (Z - 2.0 / r) *
               std::exp(-(Z * r));
    cor_U100 = -std::pow(Z, 5.0 / 2) / std::sqrt(M_PI) * (x + y + z) / r *
               std::exp(-(Z * r));
}

// Wavefunction of isolated electron in 200 state
double U200(const double x,
            const double y,
            const double z,
            const uint32_t Z)
{
    double r;
    double psi;

    r = std::sqrt(x * x + y * y + z * z);
    psi = std::pow(Z, 3.0 / 2) / std::sqrt(32 * M_PI) * (2 - r * Z) *
          std::exp(-r * Z / 2);

    return psi;
}

// Derivation of U200 wavefunction
void dU200(double &lap_U200, double &cor_U200,
           const double x, 
           const double y,
           const double z,
           const uint32_t Z)
{
    // lap_U100 is the laplacian of U100
    // cor_U100 is the sum of the differenzation of U100
    double r;
    
    r = std::sqrt(x * x + y * y + z * z);

    lap_U200 = std::pow(Z, 5.0 / 2) / std::sqrt(32 * M_PI) *
               std::exp(-(Z * r) / 2) *
               ((1 + Z * (1 - r / 2)) * (Z / 2 - 2.0 / r) + Z / 2.0);
    cor_U200 = -std::pow(Z, 5.0 / 2) / std::sqrt(32 * M_PI) * (x + y + z) / r *
               (1 + Z * (1 - r / 2)) * std::exp(-(Z * r) / 2);
}

void Jastrow(double &f_j_d, double &f_j_dx,
             double &f_j_dy, double &f_j_dz, double &f_j_l,
             const double x_ij, const double y_ij,
             const double z_ij, const double r_ij,
             const double beta)
{
    // Sum of first derivation of Jastrow function
    f_j_d += (x_ij + y_ij + z_ij) / (r_ij + beta * r_ij) /
            (r_ij + beta * r_ij) / r_ij;

    // Derivation of Jastrow in x,y and z
    f_j_dx += x_ij / (1 + beta * r_ij) / (1 + beta * r_ij) / r_ij;
    f_j_dy += y_ij / (1 + beta * r_ij) / (1 + beta * r_ij) / r_ij;
    f_j_dz += z_ij / (1 + beta * r_ij) / (1 + beta * r_ij) / r_ij;

    // Laplacian of Jastrow function
    f_j_l += 2.0 / r_ij / std::pow(1 + beta * r_ij, 3.0);
}

// Probability density function
void PSI2_tot(double &psi, double &psi_b,
              double &psi_1, double &psi_1_b,
              const std::vector<double> *x,
              const std::vector<double> *y,
              const std::vector<double> *z,
              const std::vector<double> *x1,
              const std::vector<double> *y1,
              const std::vector<double> *z1,
              const uint32_t , const double Lb,
              const double Beta, const double Beta1)
{
    double x_ij;
    double y_ij;
    double z_ij;
    double r_ij;
    constexpr uint32_t N = 4;
    uint32_t Z{N};

    // For two different beta
    double psi1{1};
    double psi2{1};
    double exp_j{0};
    double exp_j_b{0};

    // For the different position
    double psi11{1};
    double psi12{1};
    double exp_j1{0};
    double exp_j_b1{0};
    
    // first one indicate x0(0) or x1(1)
    // second one indicate U100(0) or U200(1)
    // third indicate the particle
    double coeff[2][2][N];

    for (uint32_t i{0}; i < N; i++)
    {
        // Control if a particle is go out from the box side Lb
        if (std::fabs((*x)[i]) > Lb || std::fabs((*y)[i]) > Lb || std::fabs((*z)[i]) > Lb)
        {
            psi1 *= 0;
            psi2 *= 0;
            // std::cout << "dentro 0" << std::endl;
        }
        if (std::fabs((*x1)[i]) > Lb || std::fabs((*y1)[i]) > Lb || std::fabs((*z1)[i]) > Lb)
        {
            psi11 *= 0;
            psi12 *= 0;
            // std::cout << "dentro 1" << std::endl;
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

                exp_j1 += 1.0 / 2 * r_ij / (1 + Beta * r_ij);
                exp_j_b1 += 1.0 / 2 * r_ij / (1 + Beta1 * r_ij);
            }
        }
        coeff[0][0][i]=U100((*x)[i], (*y)[i], (*z)[i], Z);
        coeff[0][1][i]=U200((*x)[i], (*y)[i], (*z)[i], Z);

        coeff[1][0][i]=U100((*x1)[i], (*y1)[i], (*z1)[i], Z);
        coeff[1][1][i]=U200((*x1)[i], (*y1)[i], (*z1)[i], Z);
    }

    // Position x0 but different beta values
    psi1 *= 0.5 * std::exp(exp_j);
    psi2 *= 0.5 * std::exp(exp_j_b);

    // Position x1 but different beta values
    psi11 *= 0.5 * std::exp(exp_j1);
    psi12 *= 0.5 * std::exp(exp_j_b1);
    
    // For x0
    psi1 *= (coeff[0][0][0] * coeff[0][1][2] - coeff[0][0][2] * coeff[0][1][0]) *
            (coeff[0][0][1] * coeff[0][1][3] - coeff[0][0][3] * coeff[0][1][1]);
    psi2 *= (coeff[0][0][0] * coeff[0][1][2] - coeff[0][0][2] * coeff[0][1][0]) *
            (coeff[0][0][1] * coeff[0][1][3] - coeff[0][0][3] * coeff[0][1][1]);

    // For x1
    psi11 *= (coeff[1][0][0] * coeff[1][1][2] - coeff[1][0][2] * coeff[1][1][0]) *
             (coeff[1][0][1] * coeff[1][1][3] - coeff[1][0][3] * coeff[1][1][1]);
    psi12 *= (coeff[1][0][0] * coeff[1][1][2] - coeff[1][0][2] * coeff[1][1][0]) *
             (coeff[1][0][1] * coeff[1][1][3] - coeff[1][0][3] * coeff[1][1][1]);

    // Probability density function
    // Calcolous of probability density punction for 2 different value of Beta
    psi = psi1 * psi1;
    psi_b = psi2 * psi2;

    // Calcolous of probability density function for x1
    psi_1 = psi11 * psi11;
    psi_1_b = psi12 * psi12;
    // std::cout<<psi_1_b<<' '<<psi_1<<std::endl;
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
            const double dt, const uint32_t ,
            const double A_r)
{
    // i0 is the old position (with i =x,y,z)
    // i1 is the new position (with i =x,y,z)
    // Fi0 is the old force (with i =x,y,z)
    // Fi1 is the new force (with i =x,y,z)
    // P0 is the probability in x_n
    // P1 is the probability in x_n+1

    constexpr uint32_t N{4};
    // Auxiliary variabile
    double sum0{0};
    double sum1{0};


    struct VARI{
        double x,y,z;
        void print()
        {
            std::cout<< x << ' ' << y << ' ' << z << std::endl;
        }
    };
    
    VARI vari[2][N];
    
    double Nu;
    double term;

    // Probability variabile
    double P;
    double A;
    

    for (uint32_t i{0}; i < N; i++)
    {
        vari[0][i].x = (*x1)[i] - (*x0)[i] - dt * (*Fx0)[i];
        vari[1][i].x = (*x0)[i] - (*x1)[i] - dt * (*Fx1)[i];

        vari[0][i].y = (*y1)[i] - (*y0)[i] - dt * (*Fy0)[i];
        vari[1][i].y = (*y0)[i] - (*y1)[i] - dt * (*Fy1)[i];

        vari[0][i].z = (*z1)[i] - (*z0)[i] - dt * (*Fz0)[i];
        vari[1][i].z = (*z0)[i] - (*z1)[i] - dt * (*Fz1)[i];

        sum0 += vari[0][i].x + vari[0][i].y + vari[0][i].z;
        sum1 += vari[1][i].x + vari[1][i].y + vari[1][i].z;
    }
    // vari[0][0].print();

    // Evaluation of exponent of transition matrix
    term = std::pow(sum1, 2.0) - std::pow(sum0, 2.0);
    Nu = term / (4 * dt * A_r);

    // Estiamtion of acceptance rate
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
               const uint32_t N, const double dt,
               const double A_r)
{

    double a, b, c, d, A, B, C;
    for (uint32_t j{0}; j < N; j++)
    {
        // Box-Muller equation for the random force
        a = (double)rand() / ((double)RAND_MAX + 1);
        b = (double)rand() / ((double)RAND_MAX + 1);
        c = (double)rand() / ((double)RAND_MAX + 1);
        d = (double)rand() / ((double)RAND_MAX + 1);

        // A,B,C have a normal distrubution with variance 2*A*dt
        A = std::sqrt(-4 * A_r * dt * std::log(1 - a)) * std::cos(2 * M_PI * b);
        B = std::sqrt(-4 * A_r * dt * std::log(1 - a)) * std::sin(2 * M_PI * b);
        C = std::sqrt(-4 * A_r * dt * std::log(1 - c)) * std::cos(2 * M_PI * d);

        // Evolution of the sistem
        (*x1)[j] = (*x0)[j] + dt * (*Fx0)[j] + A;
        (*y1)[j] = (*y0)[j] + dt * (*Fy0)[j] + B;
        (*z1)[j] = (*z0)[j] + dt * (*Fz0)[j] + C;
    }
}

// Evaluate the local energy of the sistem
void LocalEnergy(double &LE, double &LE1,
                 const std::vector<double> *x,
                 const std::vector<double> *y,
                 const std::vector<double> *z,
                 const uint32_t , const double beta,
                 const double beta1)
{
    double LocalE{0};
    double LocalE_beta{0};
    double r;
    double xx, yy, zz;


    constexpr uint32_t N=4;

    struct F_J{
        double d,dx,dy,dz,l;
    };

    F_J f_j_s[N];
    F_J f_j_s1[N];
    

    struct U_der{
        double U_lap1,U_cor1,U_lap2,U_cor2;
        void print(){
            std::cout << U_lap1 << ' ' << U_cor1 << ' ' << U_lap2 << ' ' << U_cor2;
        }
    };

    U_der u_der[N];

    struct U{
        double U100,U200;
    };

    U u[N];

    double T{0};
    double T_beta{0};

    double T0{0};
    double T0_beta{0};

    double T1{0};
    double T1_beta{0};

    double T2{0};
    double T2_beta{0};

    double T3{0};
    double T3_beta{0};

    double D0{1};
    double D1{1};

    for (uint32_t i{0}; i < N; i++)
    {
        f_j_s[i].d = 0;
        f_j_s1[i].d = 0;

        f_j_s[i].dx = 0;
        f_j_s1[i].dx = 0;

        f_j_s[i].dy = 0;
        f_j_s1[i].dy = 0;

        f_j_s[i].dz = 0;
        f_j_s1[i].dz = 0;

        f_j_s[i].l = 0;
        f_j_s1[i].l = 0;


        // Coulomb term
        r = std::sqrt((*x)[i] * (*x)[i] + (*y)[i] * (*y)[i] + (*z)[i] * (*z)[i]);

        LocalE += -(N / r);

        for (uint32_t j{0}; j < N; j++)
        {
            if (j != i)
            {
                //Interacation term electron-electron
                xx = (*x)[i] - (*x)[j];
                yy = (*y)[i] - (*y)[j];
                zz = (*z)[i] - (*z)[j];

                r = std::sqrt(xx * xx + yy * yy + zz * zz);

                // Elecron-electron interaction term 
                LocalE += 1.0 / 2 / r;
            
                // Calcolus of derivate Jastrow factor's for beta
                Jastrow(f_j_s[i].d, f_j_s[i].dx, f_j_s[i].dy, f_j_s[i].dz, f_j_s[i].l, xx, yy, zz, r, beta);
                // Calcolus of derivate Jastrow factor's for beta
                Jastrow(f_j_s1[i].d, f_j_s1[i].dx, f_j_s1[i].dy, f_j_s1[i].dz, f_j_s1[i].l, xx, yy, zz, r, beta1);
            }
        }
        dU100(u_der[i].U_lap1, u_der[i].U_cor1, (*x)[i], (*y)[i], (*z)[i], N);
        dU200(u_der[i].U_lap2, u_der[i].U_cor2, (*x)[i], (*y)[i], (*z)[i], N);
        u[i].U100 = U100((*x)[i], (*y)[i], (*z)[i], N);
        u[i].U200 = U200((*x)[i], (*y)[i], (*z)[i], N);

        T += f_j_s[i].l + f_j_s[i].dx * f_j_s[i].dx + f_j_s[i].dy * f_j_s[i].dy + f_j_s[i].dz * f_j_s[i].dz;
        T_beta += f_j_s1[i].l + f_j_s1[i].dx * f_j_s1[i].dx + f_j_s1[i].dy * f_j_s1[i].dy + f_j_s1[i].dz * f_j_s1[i].dz;
    }
    
    // Denominator of local energy
    D0 = u[0].U100 * u[2].U200 - u[0].U200 * u[2].U100;
    D1 = u[1].U100 * u[3].U200 - u[1].U200 * u[3].U100;
    
    // beta
    // Term of kinetic energy of particle i=1, in my paper
    T0 = (u_der[0].U_lap1 * u[2].U200 - u_der[0].U_lap2 * u[2].U100 + 2 * (u_der[0].U_cor1 * u[2].U200 - u[2].U100 * u_der[0].U_cor2) * f_j_s[0].d) / D0;
    // Term of kinetic energy of particle i=2, in my paper
    T1 = (u_der[1].U_lap1 * u[3].U200 - u_der[1].U_lap2 * u[3].U100 + 2 * (u_der[1].U_cor1 * u[3].U200 - u[3].U100 * u_der[1].U_cor2) * f_j_s[1].d) / D1;
    // Term of kinetic energy of particle i=3, in my paper
    T2 = (u_der[2].U_lap2 * u[0].U100 - u_der[2].U_lap1 * u[0].U200 + 2 * (u_der[2].U_cor2 * u[0].U100 - u[0].U200 * u_der[2].U_cor1) * f_j_s[2].d) / D0;
    // Term of kinetic energy of particle i=4, in my paper
    T3 = (u_der[3].U_lap2 * u[1].U100 - u_der[3].U_lap1 * u[1].U200 + 2 * (u_der[3].U_cor2 * u[1].U100 - u[1].U200 * u_der[3].U_cor1) * f_j_s[3].d) / D1;

    // beta1
    // Term of kinetic energy of particle i=1, in my paper
    T0_beta = (u_der[0].U_lap1 * u[2].U200 - u_der[0].U_lap2 * u[2].U100 + 2 * (u_der[0].U_cor1 * u[2].U200 - u[2].U100 * u_der[0].U_cor2) * f_j_s1[0].d) / D0;
    // Term of kinetic energy of particle i=2, in my paper
    T1_beta = (u_der[1].U_lap1 * u[3].U200 - u_der[1].U_lap2 * u[3].U100 + 2 * (u_der[1].U_cor1 * u[3].U200 - u[3].U100 * u_der[1].U_cor2) * f_j_s1[1].d) / D1;
    // Term of kinetic energy of particle i=3, in my paper
    T2_beta = (u_der[2].U_lap2 * u[0].U100 - u_der[2].U_lap1 * u[0].U200 + 2 * (u_der[2].U_cor2 * u[0].U100 - u[0].U200 * u_der[2].U_cor1) * f_j_s1[2].d) / D0;
    // Term of kinetic energy of particle i=4, in my paper
    T3_beta = (u_der[3].U_lap2 * u[1].U100 - u_der[3].U_lap1 * u[1].U200 + 2 * (u_der[3].U_cor2 * u[1].U100 - u[1].U200 * u_der[3].U_cor1) * f_j_s1[3].d) / D1;

    // Local energy 
    LocalE_beta = LocalE - 0.5 * (T0_beta + T1_beta + T2_beta + T3_beta + T_beta);
    LocalE += -0.5 * (T0 + T1 + T2 + T3 + T);

    LE = LocalE;
    LE1 = LocalE_beta;
}   