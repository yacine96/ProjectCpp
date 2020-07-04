#pragma once
#include "Heston_Model.h"
#include "Simulator.h"

// Functions for the simulation of ln X 
double variance_X_next(double V, double V_next, double delta, Heston_Model h);
double X_next_simulate(double X, double V, double V_next, double delta, Heston_Model h);
double expectation_X_next(double X, double V, double V_next, double delta, Heston_Model h);

// Simulator for ln X
class X_simulator : public Simulator {
public:
	X_simulator(int N, double T, double x0, Simulator* Var_Simulator, Heston_Model h = Heston_Model());
	virtual RandomProcess simulate() const;

private:
	int N;
	double T;
	double delta;
	double x0;
	Heston_Model h;
	Simulator* Var_Simulator;
};