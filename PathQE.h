#pragma once
#include "Heston_Model.h"
#include "Simulator.h"


// functions for the QE scheme 
double scheme_1(double phi, double m, double u);
double scheme_2(double phi, double m, double u);
double variance_simulation(double phi_c, double var, double delta, Heston_Model h);

// simulator for the QE scheme
const double phi_c0 = 1.5;  // default paramater phi_c 

class QE_simulator : public Simulator {
public:
	QE_simulator(int N, double T,  double v0, double phi_c = phi_c0, Heston_Model h = Heston_Model());
	virtual RandomProcess simulate() const; 
	double x0;

private:	
	double phi_c;
	double delta;
	int N;
	double T;
	Heston_Model h;
};