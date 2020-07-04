#pragma once
#include "Function.h"


const int default_max_iter = 100;
const double tol = 1e-6;

double newton_method(double x0, Function<double> const &f, int max_iter = default_max_iter);
double test_newton_method();


