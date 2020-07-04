#include "RandomProcess.h"
#include <math.h>
#include "X_simulation.h"


double expectation_X_next(double X, double V, double V_next, double delta, Heston_Model h)
{
       double integral_V = delta * (V + V_next) /2;
       return log(X) + (h.rho / h.sigma) * (V_next - V - h.kappa * h.theta * delta) + (h.kappa * h.rho / h.sigma - 0.5) * integral_V;
}

double variance_X_next(double V, double V_next, double delta, Heston_Model h)
{
    double integral_V = delta * (V + V_next) / 2;
    return (1 - h.rho * h.rho) * integral_V;
}

double X_next_simulate(double X, double V, double V_next, double delta, Heston_Model h)
{
    double expectation = expectation_X_next(X, V, V_next, delta, h);
    double variance = variance_X_next(V, V_next, delta, h);
    std::cout << variance << std::endl;
    double random_number = gaussian(expectation, sqrt(variance));
    return exp(random_number);
}

X_simulator::X_simulator(int N, double T, double x0, Simulator* Var_Simulator, Heston_Model h): N(N), T(T), delta(T / N), x0(x0), Var_Simulator(Var_Simulator), h(h)
{
}

RandomProcess X_simulator::simulate() const
{   
    double x = x0;
    double v = Var_Simulator->x0;
    double v_next;
    vector<double> data(N);

    // Simulation of the variance random process
    RandomProcess Var = Var_Simulator->simulate();
    
    // Simulation of X
    for (int i = 0; i < N; i++) {
        v_next = Var[i];
        x = X_next_simulate(x, v, v_next, delta, h);
        v = v_next;
        data[i] = x;
    }

    return RandomProcess(N, T, data);
}
