#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

class Complex
{	
	public:

		Complex();
		Complex(double r, double i);
		Complex(Complex const& a);
		double get_re() const;
		double get_im() const;
		double norm() const;
		double norm2() const;
		Complex conj() const;
		void operator+=(Complex const& a);
		void operator-=(Complex const& a);
		void operator*=(Complex const& a);
		void operator/=(double a);
		void operator/=(Complex const& a);
		Complex operator-() const;
		~Complex();


	private:

		double re;
		double im;
};

Complex operator+(Complex const& a, Complex const& b);
Complex operator+(Complex const& a, double b);
Complex operator+(double b, Complex const& a);
Complex operator-(Complex const& a, Complex const& b);
Complex operator-(Complex const& a, double b);
Complex operator-(double b, Complex const& a);
Complex operator*(Complex const& a, Complex const& b);
Complex operator*(Complex const& a, double b);
Complex operator*(double b, Complex const& a);
Complex operator/(Complex const& a, Complex const& b);
Complex operator/(Complex const& a, double b);
Complex operator/(double b, Complex const& a);
Complex exp(Complex const& a);
Complex sqrt(Complex const& a);
Complex ln(Complex const& a);
Complex pow(Complex const& a, double x);

// nombre j
const Complex j(0., 1.);

