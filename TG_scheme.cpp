#include "TG_scheme.h"

double fun_m(double var, double delta, Heston_Model h)
{
	return h.theta + (var - h.theta) * exp(-h.kappa * delta);
}

double fun_s2(double var, double delta, Heston_Model h)
{
	return (1 / h.kappa) * var * h.sigma2 * exp(-h.kappa * delta) * (1 - exp(-h.kappa * delta))
		+ h.kappa / 2 * h.theta * h.sigma2 * pow((1 - exp(-h.kappa * delta)), 2);
}

double fun_r(double psi)
{	
	class fun : public Function<double> {
	public:
		fun(double psi): psi(psi){}
		virtual double operator() (double r) const {
			double phi_r = phi(r);
			double Phi_r = normal_cdf(r);
			return r * phi_r + Phi_r * (1 + r * r) - (1 + psi) * pow((phi_r + r * Phi_r), 2);
		}
	private:
		double psi;
	};

	fun f(psi);
	double r0 = 0;
	double r = newton_method(r0, f);

	return r;
}

double tg_next_var(double var, double delta, Heston_Model h)
{
	double m = fun_m(var, delta, h);
	double s2 = fun_s2(var, delta, h);
	double s = sqrt(s2);
	double psi = s2 / (m * m);
	double r = fun_r(psi);
	double phi_r = phi(r);
	double Phi_r = normal_cdf(r);
	double mu = m * r / (phi_r + r * Phi_r);
	double sigma = s / sqrt(psi) / (phi + r * Phi_r);
	double Z = unif(0, 1);

	double res = mu + sigma * Z;

	if (res > 0)
		return res;
	
	return 0.0;
}
