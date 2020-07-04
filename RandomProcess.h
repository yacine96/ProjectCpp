#pragma once
#include<cstdlib>
#define _USE_MATH_DEFINES
#include<math.h>
#include<iostream>
#include<vector>
using namespace std;


class RandomProcess
{
public:
	RandomProcess();
	RandomProcess(int N, double T);
	RandomProcess(RandomProcess const& X);
	RandomProcess(int N, double T, vector<double> v);
	virtual ~RandomProcess();
	void operator+=(double x);
	void operator+=(RandomProcess const& X);
	void operator*=(double x);
	void operator*=(RandomProcess const& X);
	void operator-=(double x);
	void operator-=(RandomProcess const& X);
	void operator/=(double x);
	void operator/=(RandomProcess const& X);
	RandomProcess operator-() const; 
	RandomProcess inv() const;
	double operator[](int i) const;
	double operator()(double t) const; 
	void display() const;
	void exponential();
	void init_mb();
	RandomProcess diff() const;
	RandomProcess diff(int k) const;
	RandomProcess cumsum() const;
	
	
protected:
	int N;
	double T;
	vector<double> data;
};


// Operations on random proccesses
RandomProcess operator+ (RandomProcess const& X, RandomProcess const& Y);
RandomProcess operator+ (RandomProcess const& X, double y);
RandomProcess operator+ (double y, RandomProcess const& X);
RandomProcess operator- (RandomProcess const& X, RandomProcess const& Y);
RandomProcess operator- (RandomProcess const& X, double y);
RandomProcess operator- (double y, RandomProcess const& X);
RandomProcess operator* (RandomProcess const& X, RandomProcess const& Y);
RandomProcess operator* (RandomProcess const& X, double y);
RandomProcess operator* (double y, RandomProcess const& X);
RandomProcess operator/ (RandomProcess const& X, RandomProcess const& Y);
RandomProcess operator/ (RandomProcess const& X, double y);
RandomProcess operator/ (double y, RandomProcess const& X);
RandomProcess exp(RandomProcess const& X);

// Functions to simulate gaussian law and uniform law
double unif(double min, double max);
double unif();
double gaussian(double mu, double sigma);
double gaussian();

// cumulative distribution function for standardized centered normal law  
double normal_cdf(double x);

// inverse cumulative distribution fonction for standardized centered normal law 
double inverse_normal_cdf(double quantile);

// Simulate Brownian Motion
RandomProcess brownian_motion(int N, double T);

// Create a Brownian Motion Y correlated at rho% with Brownian Motion X 
RandomProcess correlated_BM(RandomProcess const& X, double rho);

// Test function
void test_RandomProcess();