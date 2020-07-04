#include "Heston_Model.h"
#include "Heston_Model.h"

Heston_Model::Heston_Model() : mu(mu0), kappa(kappa0), theta(theta0), sigma(sigma0), sigma2(sigma20), rho(rho0), r(r0)
{
}

Heston_Model::Heston_Model(double mu, double kappa, double theta, double sigma, double rho, double r) :
	mu(mu), kappa(kappa), theta(theta), sigma(sigma), sigma2(sigma* sigma), rho(rho), r(r)
{
}

Heston_Model::~Heston_Model()
{
}
