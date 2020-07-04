#pragma once
class Heston_Model
{
public:
	Heston_Model();
	Heston_Model(double mu, double kappa, double theta, double sigma, double rho, double r);
	~Heston_Model();
	double mu;
	double kappa;
	double theta;
	double sigma;
	double sigma2;
	double rho;
	double r;
};

// default parameters of heston model
const double mu0 = .03;
const double kappa0 = .50;
const double theta0 = .04;
const double sigma0 = 1.;
const double sigma20 = sigma0 * sigma0;
const double rho0 = .25;
const double r0 = 0.01;

