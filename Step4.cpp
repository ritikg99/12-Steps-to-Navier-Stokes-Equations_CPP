#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Defining the function phi
double phi(double x, double t, double nu) 
{  
    double term1 = std::exp(-std::pow(x-4 * t, 2)) / (4*nu*(t+1));
    double term2 = std::exp(-std::pow(x - 4*t - 2*3.141592654, 2) / (4*nu*(t+1)));
    return term1 + term2;
}

// // Function for computing u
// std::vector<double> initial_condition(const std::vector<double>& x, double t, double nu, double dx)
// { 
//     // std::vector<double> derivative(x.size()); 

//     // for (int i=1; i<x.size(); i++){
//     //     derivative[i] = (phi(x[i+1], t, nu) + phi(x[i-1], t, nu)) / (2*dx);
//     // }

//     // // Calcualting derivative at boundaries
//     // derivative[0] = (phi(x[1], t, nu) - phi(x[0], t, nu)) / dx;
//     // derivative[x.size() - 1] = (phi(x[x.size() - 1], t, nu) - phi(x[x.size() - 2], t, nu)) / dx;

//     std::vector<double> dphi(x.size()); 
//     for (int i=0; i<x.size(); i++){
//         dphi[i] = (exp(-pow(2*3.141592 + 4*t - x[i], 2) / (4*nu*(1 + t))) * (4*3.141592 + 8*t - 2*x[i])) / (4*nu*(1 + t)) + (exp(-pow(4*t - x[i], 2) / (4*nu*(1 + t))) * (8*t - 2*x[i])) / (4*nu*(1 + t));
//     }
//     //dphi = (exp(-pow(2*3.141592 + 4*t - x, 2) / (4*v*(1 + t))) * (4*3.141592 + 8*t - 2*x[i])) / (4*v*(1 + t)) + (exp(-pow(4*t - x, 2) / (4*v*(1 + t))) * (8*t - 2*x)) / (4*v*(1 + t));

//     // derivative[0] = -32;
//     // derivative[x.size() - 1] = -32;

//     std::vector<double> u(x.size());
//     for (int i=0; i<x.size(); i++){
//         u[i] = 4 - (2 * nu * dphi[i]/(phi(x[i], t, nu) * dx)) ;
//     }
//     plt::plot(x, u);

//     return u;
// }

// Function for computing u
double initial_condition(double t, double v, double x)
{ 
    double phi = exp(-pow(4*t - x, 2) / (4*v*(1 + t))) + exp(-pow(2*M_PI + 4*t - x, 2) / (4*v*(1 + t)));

    double dphi = (exp(-pow(2*M_PI + 4*t - x, 2) / (4*v*(1 + t))) * (4*M_PI + 8*t - 2*x)) / (4*v*(1 + t)) + (exp(-pow(4*t - x, 2) / (4*v*(1 + t))) * (8*t - 2*x)) / (4*v*(1 + t));

    double u = 4 - (2 * dphi * v) / phi;

    return u;
}

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    if (num <= 0) {
        return result; // Return an empty vector if num is non-positive
    }
    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

int main(int argc, char* argv[])
{
    int nx = 101;    // Number of dots
    double dx = 2.0*M_PI/(nx - 1);   // Local spacing
    int nt = 100;   // Number of timesteps
    double nu = 0.07;
    
    double dt = dx * nu;

    std::vector<double> un(nx);

    // Initialize the matrix with ones
    std::vector<double> u(nx,0.0);

    std::vector<double> x(nx);
    for (int i = 0; i<nx; i++)
    {
        x[i] = 2*M_PI*i/ (nx-1);
        u[i] = initial_condition(0, nu, x[i]);
    }

    double start = 0.0;
    double end = 2 * M_PI; // Use M_PI from <cmath>

    // std::vector<double> u = initial_condition(x, 0.0, nu, dx);

    for (int n=0; n<nt; n++){
        // double un[nx];
        for (int i = 0; i < nx; ++i){
            un[i] = u[i];
            u[i] = initial_condition(0, nu, x[i]);
        }

        for (int i=1; i<nx-1; i++) {
            u[i] = un[i] - (un[i] * dt/dx * (un[i]-un[i-1])) + (nu * dt/(dx*dx) * (un[i+1] - 2*un[i] + un[i-1])); 

            u[0] = un[0] - un[0] * dt/dx * (un[i] - un[i-1]) + nu * dt/(dx*dx) * (un[1] - 2 * un[0] + un[-2]);
            u[nx - 1] = u[0];
        }

        

        // if (n==0){
            plt::plot(x, u);
            plt::xlabel("x");
            plt::ylabel("u");
            // plt::xlim(-1.0,15.0);
            // plt::ylim(-40.0,10.0);
            plt::pause(0.03);
            plt::clf();
        // }
        
    }

    plt::show();  
    return 0;
}

