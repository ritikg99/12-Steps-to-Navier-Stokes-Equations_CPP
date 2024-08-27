#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(int argc, char* argv[])
{
    int nx = 41;    // Number of dots
    double dx = 2.0/(nx - 1);   // Local spacing
    int nt = 40;   // Number of timesteps
    double nu = 0.3;
    double sigma = 0.2;  // Courant number (the new dt will be calculated based on this)
    
    double dx2 = dx*dx;
    double dt = sigma * dx*dx / nu;


    std::vector<double> un(nx);

    // Initialize the matrix with ones
    std::vector<double> u(nx,0.0);

    std::vector<double> x(nx);
    for (int i = 0; i<nx; i++)
    {
        x[i] = i * dx;
    }
    
    // for (int i=9; i<18; i++)
    // {
    //     u[i] = 1.0;
    // }
    for (int i = 0; i < nx; i++) {
        x[i] = (5.0 * i) / (nx - 1);
        u[i] = (x[i] >= 0.5 && x[i] <= 1) ? 1 : 0;
    }
    plt::plot(x,u);


    for (int n=0; n<nt; n++){
        // double un[nx];
        for (int i = 0; i < nx; ++i){
            un[i] = u[i];
        }

        //u[0] = un[0];
        double dx2 = dx*dx;

        for (int i=0; i<nx; i++) {
            u[i] = un[i] + nu * dt/dx2 * (un[i+1] - 2*un[i] + un[i-1]);       
        }

        plt::plot(x, u);
        plt::xlabel("x");
        plt::ylabel("u");
        plt::xlim(-1.0,3.0);
        plt::ylim(-0.1,2.1);
        plt::pause(0.1);
        plt::clf();
    }

    plt::show();  
}
