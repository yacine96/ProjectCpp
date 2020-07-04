#pragma once
#include "RandomProcess.h"
#include "Function.h"

class Simulator
{
public:
	virtual RandomProcess simulate() const = 0;
	double monte_carlo_computation(int n_simulation, Function<double> const& f = Identity());
	double x0;

protected:
	int N;
};

