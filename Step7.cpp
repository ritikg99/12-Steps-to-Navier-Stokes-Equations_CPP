#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(int argc, char* argv[])
{
    int nx = 81;    // Number of dots in x
    int ny = 81;    // Number of dots in y
    double dx = 2.0/(nx - 1);   // Local spacing in x
    double dy = 2.0/(ny - 1);   // Local spacing in y
    int nt = 10;   // Number of timesteps
    double sigma = 0.25;
    int c = 1;  // Velocity
    int nu = .05; 
    double dt = sigma*dx*dy/ nu; // Time step size


    std::vector<double> x(nx);
    std::vector<double> y(ny);

    for (int i = 0; i<nx; i++)
    {
        x[i] = i * dx;
    }
        
    for (int i = 0; i<ny; i++)
    {
        y[i] = i * dy;
    }

    // Initialize the matrix with ones
    std::vector<std::vector<double>> u(nx, std::vector<double> (ny,1.0));

    // Matrix for u at n time step
    std::vector<std::vector<double>> un(nx, std::vector<double>(ny));


    for (int i=19; i<(nx/2); i++)
    {
        for (int j=19; j<(ny/2); j++){
            u[i][j] = 2.0;
        }
    }

    plt::ion();

    // Looping for nt time steps
    for (int n=0; n<nt; n++){
        
        // Update the un matrix to the 
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j < ny; j++){
                un[i][j] = u[i][j];
            }
        }

        for (int i=1; i<nx-1; i++) {
            for (int j = 1; j<ny-1; j++){
                u[i][j] = un[i][j] + nu*dt/(dx*dx)*(un[i+1][j] - 2*un[i][j] + un[i][j-1]) + nu*dt/(dy*dy)*(un[i][j+1] - 2*un[i][j] + un[i][j-1]);
            }
        }

        // Boundary conditions
        for(int i=0; i<nx; i++) u[0][i] = 1.0;
        for(int i=0; i<ny; i++) u[i][0] = 1.0;

        for (int i = 0; i < nx; i++) u[i][ny-1] = 1.0;
        for (int i = 0; i < ny; i++) u[nx-1][i] = 1.0;

        std::vector<std::vector<double>> X(nx, std::vector<double>(ny));
        std::vector<std::vector<double>> Y(nx, std::vector<double>(ny));
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                X[i][j] = i * dx;
                Y[i][j] = j * dy;
        }
    }

        // , {}, animatedFig const long animatedFig = plt::figure(1);
        
        plt::clf();
        // plt::surface_plot(X, Y, u);
        // plt::plot(X,u)
        plt::xlabel("x");
        plt::ylabel("y");
        plt::pause(0.1);
    }

    plt::show();
    return 0;
}
