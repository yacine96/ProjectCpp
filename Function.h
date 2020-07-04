#pragma once
#include <math.h>

template <typename T>
class Function
{
	public:
		
		virtual T operator()(double w) const = 0;
		T diff1(double w) const;
		T diff2(double w) const;
		virtual ~Function();

	private:
		double eps = 1e-4;
};

template <typename T>
T Function<T>::diff1(double w) const
{
	T a = this->operator()(w + eps);
	T b = this->operator()(w);
	return (a - b) / eps;
}

template <typename T>
T Function<T>::diff2(double w) const
{
	T a = this->operator()(w + eps);
	T b = this->operator()(w + eps / 2);
	T c = this->operator()(w);

	return ((a + c) - 2 * b) / (eps * eps / 4);
}

template <typename T>
Function<T>::~Function()
{
}

class Identity : public Function<double> {
public:
	virtual double operator() (double x) const {
		return x;
	};
};