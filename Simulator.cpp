#include "Simulator.h"

double Simulator::monte_carlo_computation(int n_simulation, Function<double> const& f)
{	
	double res = 0;
	for (int i = 0; i < n_simulation; i++)
		res += f(this->simulate()[N - 1]);

	return res / n_simulation;
}
