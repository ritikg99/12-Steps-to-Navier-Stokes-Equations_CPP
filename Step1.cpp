#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(int argc, char* argv[])
{
    int nx = 41;    // Number of dots
    double dx = 2.0/(nx - 1);   // Local spacing
    int nt = 50;   // Number of timesteps
    double dt = 0.025;  // Time step size
    int c = 1;  // Velocity


    std::vector<double> rows(nx, 1.0), un(nx);

    // Initialize the matrix with ones
    std::vector<double> u(nx,1.0);

    std::vector<double> x(nx);
    for (int i = 0; i<nx; i++)
    {
        x[i] = i * dx;
    }
    
    for (int i=9; i<(nx/2); i++)
    {
        u[i] = 4.0;
    }


    for (int n=0; n<nt; n++){
        // double un[nx];
        for (int i = 0; i < nx; ++i){
            un[i] = u[i];
        }

        for (int i=0; i<nx; i++) {
            u[i] = un[i] - c*dt/dx*(un[i]-u[i-1]);
        }

        plt::plot(x, u);
        plt::xlabel("x");
        plt::ylabel("y");
        plt::pause(0.1);

    }

    plt::show();
    
}
