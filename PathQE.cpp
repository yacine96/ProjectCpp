#include <math.h>
#include "PathQE.h"
#include "RandomProcess.h"


double scheme_1(double phi, double m, double u)
{
    double b2 = 2 / phi - 1 + sqrt(2 / phi) * sqrt(2 / phi -1);
    double a = m / (1 + b2);
    double Z = inverse_normal_cdf(u);

    return a * pow((Z + sqrt(b2)), 2);

}

double scheme_2(double phi, double m, double u)
{
    double p = (phi - 1)/(phi + 1);
    double beta = (1 - p) / m;
    double inv_distrib = 0;
    if (u > p) {
        inv_distrib = (1 / beta) * log((1 - p)/(1 - u));
    }
    return inv_distrib;
}

double variance_simulation(double phi_c, double var, double delta, Heston_Model h)
{
    double next_var;
    double m = h.theta + (var - h.theta) * exp(-h.kappa * delta);
    double s2 = (1 / h.kappa) * var * h.sigma2 * exp(-h.kappa * delta) * (1 - exp(-h.kappa * delta))
                + h.kappa / 2 * h.theta * h.sigma2 * pow((1 - exp(-h.kappa * delta)), 2);

    double phi = s2 / pow(m,2);
    double random_number = unif(0, 1);
    if (phi_c < phi)
    {
        next_var = scheme_2(phi, m, random_number);
    }
    else
    {
        next_var = scheme_1(phi, m, random_number);
    }
    return next_var;
}

QE_simulator::QE_simulator(int N, double T, double v0, double phi_c, Heston_Model h): x0(v0), phi_c(phi_c), delta(T / N), N(N), T(T), h(h)
{
}

RandomProcess QE_simulator::simulate() const
{   
    double v = x0;
    vector<double> data(N);
    for (int i = 0; i < N; i++) {
        v = variance_simulation(phi_c, v, delta, h);
        data[i + 1] = v;
    }
    return RandomProcess(N, N * delta, data);
}
