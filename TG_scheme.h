#pragma once
#include "Simulator.h"
#include "Newton_method.h"
#include "Heston_Model.h"

double fun_m(double var, double delta, Heston_Model h);
double fun_s2(double var, double delta, Heston_Model h);



double fun_psi(double var, double delta, Heston_Model h);
double fun_r(double psi);
double fun_mu(double r);
double fun_sigma(double r);
