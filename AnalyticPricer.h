#pragma once
#include "Function.h"
#include "Heston_Model.h"
#include "Complex.h"

class AnalyticPricer
{
	public:
		AnalyticPricer(double T, int N);
		AnalyticPricer(double T, int N, Heston_Model hm);
		double compute_swap_price(double v0);

	private:
		Heston_Model hm;
		double T;
		int N;
};

class Fun_D : public Function<Complex>
{
public:
	Fun_D(double tau, Heston_Model hm);
	virtual Complex operator()(double w) const;
	virtual ~Fun_D();

private:
	double tau;
	Heston_Model hm;
};

class Fun_C : public Function<Complex>
{
public:
	Fun_C(double tau, Heston_Model hm);
	virtual Complex operator()(double w) const;
	virtual ~Fun_C();

private:
	double tau;
	Heston_Model hm;
};


// parameters for functions C and D
Complex fun_a(double w, double kappa, double rho, double sigma);
Complex fun_delta(double w, double kappa, double rho, double sigma);
Complex fun_b(double w, double kappa, double rho, double sigma);

// functions C and D
Complex fun_d(double w, double tau, Heston_Model hm);
Complex fun_c(double w, double tau, Heston_Model hm);

// 1st & 2nd Derivatives of C & D with respect to w 
Complex fun_d1(double w, double tau, Heston_Model hm);
Complex fun_c1(double w, double tau, Heston_Model hm);
Complex fun_d2(double w, double tau, Heston_Model hm);
Complex fun_c2(double w, double tau, Heston_Model hm);

// function f + 1st & 2nd derivatives
Complex fun_f(double w, double tau, Heston_Model hm, double v, double x);
Complex fun_f1(double w, double tau, Heston_Model hm, double v, double x);
Complex fun_f2(double w, double tau, Heston_Model hm, double v, double x);

// function G 
Complex G(double T, int N, Heston_Model hm, double v);

// functions c, q_tilde, W
double fun_c_i(double t, Heston_Model hm);
double q_tilde(Heston_Model hm);
double fun_W(double t, Heston_Model hm, double v0);
Complex fun_U(double ti, double tj, Heston_Model hm, double v0);

// solution
Complex solution(double T, int N, int i, Heston_Model hm, double v0);
Complex full_solution(double T, int N, Heston_Model hm, double v0);

// closed formula when N -> +inf 
double closed_formula(double T, Heston_Model hm, double v0);