#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::vector<std::vector<double>> Laplace(std::vector<std::vector<double>> p, double tolerance, double dx, double dy, int nx, int ny){
    
    double l1norm = 1.0;
    std::vector<std::vector<double>> pn = p;
    
    while (l1norm > tolerance){
               
        pn = p;
        l1norm = 0.0;

        // Obtianing the value of p at "i"th iteration
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                p[i][j] = (dy * dy * (pn[i + 1][j] + pn[i - 1][j]) +
                           dx * dx * (pn[i][j + 1] + pn[i][j - 1])) /
                          (2 * (dx * dx + dy * dy));
                l1norm += std::fabs(p[i][j] - pn[i][j]);
            }
        }

        // Boundary conditions
        for (int i=0; i<nx; i++){
            p[0][i] = 0;
            p[nx-1][i] = 2;
            p[i][0] = p[i][1];
            p[i][ny-1] = p[i][ny-2];
        }
    }
    return p;
}



int main(int argc, char* argv[])
{
    int nx = 81;    // Number of dots in x
    int ny = 81;    // Number of dots in y
    double dx = 2.0/(nx - 1);   // Local spacing in x
    double dy = 2.0/(ny - 1);   // Local spacing in y
    int c = 1;  // Velocity

    double tolerance = .01;

    std::vector<double> x(nx);
    std::vector<double> y(ny);

    for (int i = 0; i<nx; i++){
        x[i] = i * dx;
    }
        
    for (int i = 0; i<ny; i++){
        y[i] = i * dy;
    }

    std::vector<std::vector<double>> p(nx, std::vector<double> (ny,0.0));

    std::vector<std::vector<double>> X(nx, std::vector<double>(ny));
    std::vector<std::vector<double>> Y(nx, std::vector<double>(ny));
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            X[i][j] = i * dx;
            Y[i][j] = j * dy;
        }
    }

    // Boundary conditions
    for (int i=0; i<nx; i++){
        p[0][i] = 0;
        p[nx-1][i] = 2;
        p[i][0] = p[i][1];
        p[i][ny-1] = p[i][ny-2];
    }
        
    // Calling the Laplace function
    p = Laplace(p, tolerance, dx, dy, nx, ny);

    plt::clf();
    plt::xlabel("x");
    plt::ylabel("y");
    plt::plot_surface(X,Y,p);
    // plt::pause(0.1);
    plt::show();

    return 0;
}