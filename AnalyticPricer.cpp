#include "AnalyticPricer.h"
#include "AnalyticPricer.h"

//
AnalyticPricer::AnalyticPricer(double T, int N): T(T), N(N), hm(Heston_Model())
{
}

AnalyticPricer::AnalyticPricer(double T, int N, Heston_Model hm): T(T), N(N), hm(hm)
{
}

double AnalyticPricer::compute_swap_price(double v0)
{
	return full_solution(T, N, hm, v0).get_re();
}

Fun_D::Fun_D(double tau, Heston_Model hm) : tau(tau), hm(hm)
{
}

Complex Fun_D::operator()(double w) const
{
	double kappa = hm.kappa;
	double rho = hm.rho;
	double sigma = hm.sigma;
	double sigma2 = hm.sigma2;
	Complex a(fun_a(w, kappa, rho, sigma));
	Complex b(fun_b(w, kappa, rho, sigma));
	Complex g = (a - b) / (a + b);

	return (a - b) / sigma2 * (1 - exp(-tau * b)) / (1 - g * exp(-tau * b));
}

Fun_D::~Fun_D()
{
}

Fun_C::Fun_C(double tau, Heston_Model hm) : tau(tau), hm(hm)
{
}

Complex Fun_C::operator()(double w) const
{
	double kappa = hm.kappa;
	double rho = hm.rho;
	double sigma = hm.sigma;
	double sigma2 = hm.sigma2;
	double mu = hm.mu;
	double theta = hm.theta;
	double r = hm.r;
	Complex a = fun_a(w, kappa, rho, sigma);
	Complex b = fun_b(w, kappa, rho, sigma);
	Complex g = (a - b) / (a + b);

	return (mu * j * w - r) * tau + kappa * theta / sigma2 * ((a - b) * tau - 2 * ln((1 - g * exp(-tau * b)) / (1 - g)));
}

Fun_C::~Fun_C()
{
}


Complex fun_a(double w, double kappa, double rho, double sigma)
{
	return kappa - rho * sigma * j * w;
}

Complex fun_delta(double w, double kappa, double rho, double sigma)
{
	Complex a(fun_a(w, kappa, rho, sigma));
	double sigma2 = sigma * sigma;

	return a * a + sigma2 * (j * w + w * w);
}

Complex fun_b(double w, double kappa, double rho, double sigma)
{
	return sqrt(fun_delta(w, kappa, rho, sigma));
}

Complex fun_d(double w, double tau, Heston_Model hm)
{
	return Fun_D(tau, hm)(w);
}

Complex fun_c(double w, double tau, Heston_Model hm)
{
	return Fun_C(tau, hm)(w);
}

Complex fun_d1(double w, double tau, Heston_Model hm)
{
	return Fun_D(tau, hm).diff1(w);
}

Complex fun_c1(double w, double tau, Heston_Model hm)
{
	return Fun_C(tau, hm).diff1(w);
}

Complex fun_d2(double w, double tau, Heston_Model hm)
{
	return Fun_D(tau, hm).diff2(w);
}

Complex fun_c2(double w, double tau, Heston_Model hm)
{
	return Fun_C(tau, hm).diff2(w);
}

Complex fun_f(double w, double tau, Heston_Model hm, double v, double x)
{
	Complex C = fun_c(w, tau, hm);
	Complex D = fun_d(w, tau, hm);

	return exp(C + D * v + j * w * x);
}

Complex fun_f1(double w, double tau, Heston_Model hm, double v, double x)
{
	Complex C1 = fun_c1(w, tau, hm);
	Complex D1 = fun_d1(w, tau, hm);

	return (C1 + D1 * v + j * x) * fun_f(w, tau, hm, v, x);
}

Complex fun_f2(double w, double tau, Heston_Model hm, double v, double x)
{
	Complex C1 = fun_c1(w, tau, hm);
	Complex C2 = fun_c2(w, tau, hm);
	Complex D1 = fun_d1(w, tau, hm);
	Complex D2 = fun_d2(w, tau, hm);

	return (C2 + D2 * v + (C1 + D1 * v + j * x) * (C1 + D1 * v + j * x)) * fun_f(w, tau, hm, v, x);
}

Complex G(double T, int N, Heston_Model hm, double v)
{
	double tau = T / N;
	double w = 0;
	double v2 = v * v;
	double r = hm.r;
	Complex D = fun_d(w, tau, hm);
	Complex D1 = fun_d1(w, tau, hm);
	Complex D2 = fun_d2(w, tau, hm);
	Complex C = fun_c(w, tau, hm);
	Complex C1 = fun_c1(w, tau, hm);
	Complex C2 = fun_c2(w, tau, hm);

	Complex res = -(D1 * D1 * v2 + (2 * C1 * D1 + D2) * v + C1 * C1 + C2);

	return res;
}

double fun_c_i(double t, Heston_Model hm)
{
	double kappa = hm.kappa;
	double sigma2 = hm.sigma2;
	return 2 * kappa / (sigma2 * (1 - exp(-kappa * t)));
}

double q_tilde(Heston_Model hm)
{
	double kappa = hm.kappa;
	double theta = hm.theta;
	double sigma2 = hm.sigma2;
	return 2 * kappa * theta / sigma2;
}

double fun_W(double t, Heston_Model hm, double v0)
{
	double kappa = hm.kappa;
	return fun_c_i(t, hm) * v0 * exp(-kappa * t);
}

Complex fun_U(double ti, double tj, Heston_Model hm, double v0)
{
	// tj = t_(i-1)
	Complex D1 = fun_d1(0, ti - tj, hm);
	Complex D2 = fun_d2(0, ti - tj, hm);
	Complex C1 = fun_c1(0, ti - tj, hm);
	Complex C2 = fun_c2(0, ti - tj, hm);
	double c = fun_c_i(tj, hm);
	double W = fun_W(tj, hm, v0);
	double r = hm.r;
	Complex U1 = D1 * D1 * ((q_tilde(hm) + 2 * W) + (q_tilde(hm) + W) * (q_tilde(hm) + W));
	U1 /= c * c;
	Complex U2 = (2 * C1 * D1 + D2) * (q_tilde(hm) + W) / c;
	Complex U3 = C1 * C1 + C2;

	Complex res = -(U1 + U2 + U3);

	return res;
}

Complex solution(double T, int N, int i, Heston_Model hm, double v0)
{
	double r = hm.r;

	if (i == 1) {
		double t = T / N;
		return exp(-r * t) * G(T, N, hm, v0);
	}

	else {
		double ti = T / N * i;
		double tj = T / N * (i - 1.);
		return fun_U(ti, tj, hm, v0);
	}
}

Complex full_solution(double T, int N, Heston_Model hm, double v0)
{
	Complex res(0, 0);
	for (int i = 1; i < N + 1; i++) {
		res += solution(T, N, i, hm, v0) / T;
	}
	return res;
}

double closed_formula(double T, Heston_Model hm, double v0)
{	
	double res = hm.theta + (v0 - hm.theta) * (1 - exp(-hm.kappa * T)) / (hm.kappa * T);

	return res;
}
