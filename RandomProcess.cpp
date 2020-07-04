#include "RandomProcess.h"

RandomProcess::RandomProcess() : N(1), T(1), data(0)
{	
	data = vector<double>(N);
	for (int i = 0; i < N; i++)
		data[i] = 0;
}

RandomProcess::RandomProcess(int N, double T): N(N), T(T), data(0)
{
	data = vector<double>(N);
	for (int i = 0; i < N; i++)
		data[i] = 0;
}

RandomProcess::RandomProcess(RandomProcess const& X): N(X.N), T(X.T), data(0)
{
	data = vector<double>(N);
	for (int i=0; i<X.N; i++) {
		data[i] = X.data[i];
	}
}

RandomProcess::RandomProcess(int N, double T, vector<double> v): N(N), T(T)
{
	data = vector<double>(N);
	for (int i = 0; i < N; i++)
		data[i] = v[i];
}

void RandomProcess::operator+=(double x)
{
	for (int i = 0; i < N; i++) {
		data[i] += x;
	}
}

void RandomProcess::operator+=(RandomProcess const& X)
{
	for (int i = 0; i < N; i++) {
		data[i] += X.data[i];
	}
}

void RandomProcess::operator*=(double x)
{
	for (int i = 0; i < N; i++) {
		data[i] *= x;
	}
}

void RandomProcess::operator*=(RandomProcess const& X)
{
	for (int i = 0; i < N; i++) {
		data[i] *= X.data[i];
	}
}

void RandomProcess::operator-=(double x)
{
	this->operator+=(-x);
}

void RandomProcess::operator-=(RandomProcess const& X)
{
	RandomProcess Y(X);
	Y *= -1;
	this->operator+=(Y);
}

void RandomProcess::operator/=(double x)
{
	this->operator*=(1 / x);
}

void RandomProcess::operator/=(RandomProcess const& X)
{
	for (int i = 0; i < N; i++) {
		data[i] /= X.data[i];
	}
}


RandomProcess RandomProcess::operator-() const
{	
	RandomProcess X(*this);
	X *= -1;
	return X;
}

RandomProcess RandomProcess::inv() const
{
	RandomProcess X(*this);
	for (int i = 0; i < N; i++) {
		X.data[i] = 1 / data[i];
	}

	return X;
}

double RandomProcess::operator[](int i) const
{
	return data[i];
}

double RandomProcess::operator()(double t) const
{
	int i(t / T * (double(N) - 1));
	return this->operator[](i);
}

void RandomProcess::display() const
{
	for (int i = 0; i < N; i++) {
		cout << data[i] << " ";
	}
	cout << endl;
}	

void RandomProcess::exponential()
{
	for (int i = 0; i < N; i++)
	{
		data[i] = exp(data[i]);
	}
}

void RandomProcess::init_mb()
{
	data[0] = 0;
	for (int i = 1; i < N; i++)
		data[i] = data[i-1] + gaussian(0, sqrt(T / (double(N) - 1)));
}

RandomProcess RandomProcess::diff() const
{	
	RandomProcess X(*this);
	for (int i = 0; i < N - 1; i++)
		X.data[i + 1] = data[i + 1] - data[i];
	return X;
}

RandomProcess RandomProcess::diff(int k) const
{	
	if (k == 0)
		return this->diff();
	else
	{
		RandomProcess X = this->diff();
		return X.diff(k - 1);
	}
}

RandomProcess RandomProcess::cumsum() const
{	
	RandomProcess X(*this);
	X.data[0] = data[0];
	for (int i = 1; i < N; i++)
		X.data[i] = X.data[i - 1] + data[i];
	return X;
}

RandomProcess::~RandomProcess()
{
}

RandomProcess operator+(RandomProcess const& X, RandomProcess const& Y)
{
	RandomProcess Z(X);
	Z += Y;
	return Z;
}

RandomProcess operator+(RandomProcess const& X, double y)
{
	RandomProcess Z(X);
	Z += y;
	return Z;
}

RandomProcess operator+(double y, RandomProcess const& X)
{
	return X + y;
}

RandomProcess operator-(RandomProcess const& X, RandomProcess const& Y)
{
	RandomProcess Z(X);
	Z -= Y;
	return Z;
}

RandomProcess operator-(RandomProcess const& X, double y)
{
	RandomProcess Z(X);
	Z -= y;
	return Z;
}

RandomProcess operator-(double y, RandomProcess const& X)
{
	return -(X - y);
}

RandomProcess operator*(RandomProcess const& X, RandomProcess const& Y)
{
	RandomProcess Z(X);
	Z *= Y;
	return Z;
}

RandomProcess operator*(RandomProcess const& X, double y)
{
	RandomProcess Z(X);
	Z *= y;
	return Z;
}

RandomProcess operator*(double y, RandomProcess const& X)
{
	return X * y;
}

RandomProcess operator/(RandomProcess const& X, RandomProcess const& Y)
{
	RandomProcess Z(X);
	Z /= Y;
	return Z;
}

RandomProcess operator/(RandomProcess const& X, double y)
{
	RandomProcess Z(X);
	Z /= y;
	return Z;
}

RandomProcess operator/(double y, RandomProcess const& X)
{
	return y * X.inv();
}

RandomProcess exp(RandomProcess const& X)
{	
	RandomProcess Y = RandomProcess(X);
	Y.exponential();
	return Y;
}

RandomProcess brownian_motion(int N, double T)
{
	RandomProcess X(N, T);
	X.init_mb();
	return X;
}

RandomProcess correlated_BM(RandomProcess const& X, double rho)
{
	RandomProcess Y(X);
	Y.init_mb();
	return rho * X + sqrt(1 - rho * rho) * Y;
}

double unif(double min, double max)
{
	return double(rand()) / double(RAND_MAX) * (max - min) + min;
}

double unif()
{
	return unif(0, 1);
}

double gaussian(double mu, double sigma)
{
	double u1 = unif();
	double u2 = unif();
	double g = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
	return g * sigma + mu;
}

double gaussian()
{
	return gaussian(0, 1);
}

double normal_cdf(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

	return 0.5 * (1.0 + sign * y);

}

double inverse_normal_cdf(double quantile)
{
	static double a[4] = { 2.50662823884,
					 -18.61500062529,
					  41.39119773534,
					 -25.44106049637 };

	static double b[4] = { -8.47351093090,
						  23.08336743743,
						 -21.06224101826,
						   3.13082909833 };

	static double c[9] = { 0.3374754822726147,
						0.9761690190917186,
						0.1607979714918209,
						0.0276438810333863,
						0.0038405729373609,
						0.0003951896511919,
						0.0000321767881768,
						0.0000002888167364,
						0.0000003960315187 };

	if (quantile >= 0.5 && quantile <= 0.92) {
		double num = 0.0;
		double denom = 1.0;

		for (int i = 0; i < 4; i++) {
			num += a[i] * pow((quantile - 0.5), 2 * i + 1);
			denom += b[i] * pow((quantile - 0.5), 2 * i);
		}
		return num / denom;
	}
	else if (quantile > 0.92 && quantile < 1) {
		double num = 0.0;

		for (int i = 0; i < 9; i++) {
			num += c[i] * pow((log(-log(1 - quantile))), i);
		}
		return num;
	}
	else {
		return -1.0 * inverse_normal_cdf(1 - quantile);
	}

}

void test_RandomProcess()
{
	double T = 1;
	int N = 6;
	RandomProcess X(N, T);
	cout << "Test 1" << endl;
	X.display();
	cout << "test addition" << endl;
	X += 1;
	X.display();
	cout << "test copie" << endl;
	RandomProcess Y(X);
	Y.display();
	Y += 1;
	Y.display();
	X.display();
	cout << "test addition" << endl;
	(X + Y).display();
	cout << "test multiplication" << endl;
	(X * Y).display();
	cout << "test division" << endl;
	(X / Y).display();
	cout << "test exponentielle" << endl;
	exp(X).display();
	cout << "test operator -()" << endl;
	(-X).display();
	cout << "test MB" << endl;
	RandomProcess Z = brownian_motion(N, T);
	RandomProcess Z2 = brownian_motion(N, T);
	Z.display();
	Z2.display();
	cout << "test correlated MB" << endl;
	correlated_BM(Z, .8).display();
}
