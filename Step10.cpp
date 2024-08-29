#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(int argc, char* argv[])
{
    int nx = 50;    // Number of dots in x
    int ny = 50;    // Number of dots in y
    double dx = 2.0/(nx - 1);   // Local spacing in x
    double dy = 2.0/(ny - 1);   // Local spacing in y
    int nt = 100;
    int c = 1;  // Velocity

    double tolerance = .01;

    // Initilization 
    std::vector<double> x(nx);
    std::vector<double> y(ny);

    for (int i = 0; i<nx; i++){
        x[i] = i * dx;
    }
        
    for (int i = 0; i<ny; i++){
        y[i] = i * dy;
    }

    std::vector<std::vector<double>> p(nx, std::vector<double> (ny,0.0));
    std::vector<std::vector<double>> pn = p;

    // Source terms
    std::vector<std::vector<double>> b(nx, std::vector<double> (ny,0.0));
    b[20][20] = 100;
    b[60][60] = 100;

    // Loop through time steps nt
    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++)
        {
            pn = p;
            
            p[i][j] = ((dy*dy*(pn[i+1][j] + pn[i-1][j])) + dx*dx*(pn[i][j+1] + pn[i][j-1]) 
                        - b[i][j]*dx*dx*dy*dy) / (2*(dx*dx + dy*dy));

            // Boundary conditions
            for (int i=0; i<nx; i++){
                p[0][i] = 0;
                p[ny-1][i] = 0;
                p[i][0] = 0;
                p[i][nx-1] = 0;
            }
        }
    }

    // Plotting
    std::vector<std::vector<double>> X(nx, std::vector<double>(ny));
    std::vector<std::vector<double>> Y(nx, std::vector<double>(ny));
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            X[i][j] = i * dx;
            Y[i][j] = j * dy;
        }
    }

    plt::clf();
    plt::xlabel("x");
    plt::ylabel("y");
    plt::plot_surface(X,Y,p);
    // plt::pause(0.1);
    plt::show();

    return 0;
}