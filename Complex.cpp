#include "Complex.h"
#include "Complex.h"

Complex::Complex() : re(0), im(0)
{
}

Complex::Complex(double r, double i): re(r), im(i)
{
}

Complex::Complex(Complex const& a): re(a.re), im(a.im)
{
}

double Complex::get_re() const
{
	return re;
}

double Complex::get_im() const
{
	return im;
}

double Complex::norm() const
{
	return sqrt(re * re + im * im);
}

double Complex::norm2() const
{
	return re * re + im * im;
}

Complex Complex::conj() const
{
	return Complex(re, -im);
}

void Complex::operator+=(Complex const& a)
{
	re += a.re;
	im += a.im;
}

void Complex::operator-=(Complex const& a)
{
	re -= a.re;
	im -= a.im;
}

void Complex::operator*=(Complex const& a)
{
	double temp(re);
	re = re * a.re - im * a.im;
	im = temp * a.im + im * a.re;
}

void Complex::operator/=(double a)
{
	re /= a;
	im /= a;
}

void Complex::operator/=(Complex const& a)
{
	Complex b(a.conj());
	*this /= b.norm2();
	*this *= b;
}

Complex Complex::operator-() const
{
	Complex b(*this);
	b *= Complex(-1, 0); 
	return b;
}

Complex::~Complex()
{
}

Complex operator+(Complex const& a, Complex const& b)
{	
	Complex c(a);
	c += b;
	return c;
}

Complex operator+(Complex const& a, double b)
{	
	return Complex(a) + Complex(b, 0);
}

Complex operator+(double b, Complex const& a)
{
	return a + b;
}

Complex operator-(Complex const& a, Complex const& b)
{
	Complex c(a);
	c -= b;
	return c;
}

Complex operator-(Complex const& a, double b)
{
	return Complex(a) - Complex(b, 0);
}

Complex operator-(double b, Complex const& a)
{
	return -(a - b);
}

Complex operator*(Complex const& a, Complex const& b)
{	
	Complex c(a);
	c *= b;
	return c;
}

Complex operator*(Complex const& a, double b)
{
	return Complex(a) * Complex(b, 0);
}

Complex operator*(double b, Complex const& a)
{
	return a * b;
}

Complex operator/(Complex const& a, Complex const& b)
{	
	Complex c(a);
	c /= b;
	return c;
}

Complex operator/(Complex const& a, double b)
{
	return a * (1 / b);
}

Complex operator/(double b, Complex const& a)
{
	return Complex(b, 0) / a;
}

Complex exp(Complex const& a)
{	
	double norm(exp(a.get_re()));
	Complex angle(cos(a.get_im()), sin(a.get_im()));
	return norm * angle ;
}

Complex sqrt(Complex const& a)
{	
	double im(a.get_im());
	double re(a.get_re());
	double t;

	if (re == 0) {
		t = M_PI_4;
	}

	else {
		t = atan(im / re) / 2;
	}
	return sqrt(a.norm()) * exp(Complex(0., t));
}

Complex ln(Complex const& a)
{	
	double re = log(a.norm());
	double im(0);

	if (a.get_re() == 0) {
		im = M_PI_2;
	}

	else {
		im = atan(a.get_im() / a.get_re());
	}

	return Complex(re, im);
}

Complex pow(Complex const& a, double x)
{
	return exp(ln(a) * x);
}
