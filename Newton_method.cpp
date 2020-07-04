#include "Newton_method.h"
#include <cmath>

double newton_method(double x0, Function<double> const &f, int max_iter)
{	
	int iter = 0;
	double y = f(x0);
	double x = x0;

	while ((abs(y) > tol) and (iter < max_iter)) {
		x -= f(x) / f.diff1(x);
		y = f(x);
		iter += 1;
	}

	return x;
}

double test_newton_method()
{	
	class x_square : public Function<double> {
		virtual double operator() (double x) const {
			return x * x;
		}
	};

	double x0 = 100;
	x_square fun_x_square;
	double solution = newton_method(x0, fun_x_square);

	return solution;
}
