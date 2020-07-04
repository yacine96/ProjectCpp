#include "AnalyticPricer.h"
#include <iostream>
#include "RandomProcess.h"
#include <iostream>
#include "PathQE.h"
#include "X_simulation.h"
#include "Newton_method.h"
#include <vector>

using namespace std;

int main()
{	
	//double T = 10;
	//int N = 10000;
	//Heston_Model hm;
	//AnalyticPricer Pricer(T, N, hm);
	//double v0 = theta0 * 1.5;

	//cout << Pricer.compute_swap_price(v0) << endl;
	//cout << closed_formula(T, hm, v0) << endl;

	//double x_test = X_next_simulate(1, 0, 0, 1, 0.05, 0, 0.5, 1);
	//double var_test = variance_simulation(1.5, 0.03, 0.5, 0.01, 1, 0.1);
	//cout << x_test << endl;
	//cout << var_test << endl;

	//double y = test_newton_method();
	//cout << y;
	
	//int N = 10;
	//vector<double> v(N);
	//for (int i = 0; i < N; i++)
	//	v[i] = 1;

	//cout << v[N - 1] << endl;

	return 0;
}

