#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(int argc, char* argv[])
{
    int nx = 41;    // Number of dots in x
    int ny = 41;    // Number of dots in y
    int nt = 100;   // Number of time steps
    double dx = 2.0/(nx - 1);   // Local spacing in x
    double dy = 2.0/(ny - 1);   // Local spacing in y
    int c = 1;  // Velocity

    double rho = 1.0;
    double nu = .1;
    double dt = .001;

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

    std::vector<std::vector<double>> u(nx, std::vector<double> (ny,0.0));
    std::vector<std::vector<double>> v(nx, std::vector<double> (ny,0.0));
    // std::vector<std::vector<double>> un(nx, std::vector<double> (ny,0.0));
    // std::vector<std::vector<double>> vn(nx, std::vector<double> (ny,0.0));
    std::vector<std::vector<double>> p(nx, std::vector<double> (ny,0.0));
    std::vector<std::vector<double>> b(nx, std::vector<double> (ny,0.0));

    // Lid-driven cavity condition
    for (int i = 0; i < ny; i++) {
        u[nx - 1][i] = 1.0;
    }

    for (int i=0; i<nt; i++)
    {
        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;

        // Build up b (This is the part of pressure equation in the end in brackets so that the big pressure equation is split up for better management)
        for (int i=1; i<nx-1; i++){
            for (int j=1; j<ny-1; j++)
            {
                b[i][j] = rho * (1.0/dt*(((u[i+1][j] - u[i-1][j])/(2.0*dx)) + ((v[i][j+1] - v[i][j-1])/(2.0*dy)))
                            - ((u[i+1][j] - u[i-1][j])/(2.0*dx))*((u[i+1][j] - u[i-1][j])/(2.0*dx)) 
                            - (2.0*((u[i][j+1] - u[i][j-1])/(2.0*dy))*((v[i+1][j] - v[i-1][j])/(2.0*dx))) 
                            - ((v[i][j+1] - v[i][j-1])/(2.0*dy))*((v[i][j+1] - v[i][j-1])/(2.0*dy)));
            }
        }

        std::vector<std::vector<double>> pn(nx, std::vector<double> (ny,0.0));

        // Loop through time steps nt
        for (int it=0; it<50; it++){
            pn = p;

            for (int i=1; i<nx-1; i++){
                for (int j=1; j<ny-1; j++)
                {
                    p[i][j] = (((dy*dy*(pn[i+1][j] + pn[i-1][j])) + dx*dx*(pn[i][j+1] + pn[i][j-1])) / (2.0*(dx*dx + dy*dy)))
                                - ((dx*dx*dy*dy)/(2.0*(dx*dx + dy*dy))) * (b[i][j]) ;
                }
            }

            // Boundary conditions
            // Pressure boundary conditions
            for (int i = 0; i < nx; i++) {
                p[i][0] = p[i][1];
                p[i][ny - 1] = p[i][ny - 2];
            }

            for (int j = 0; j < ny; j++) {
                p[0][j] = p[1][j];
                p[nx - 1][j] = 0.0;
            }
        }

        ///////////////// Computing u and v ///////////////// 
        for (int i=1; i<nx-1; i++){
            for (int j=1; j<ny-1; j++){
                u[i][j] = un[i][j] - un[i][j]*(dt/dx)*(un[i][j] - un[i-1][j]) - vn[i][j]*(dt/dy)*(un[i][j] - un[i][j-1])
                            - (dt/(rho*2.0*dx)) * (pn[i+1][j] - pn[i-1][j])
                            + nu*((dt/(dx*dx))*(un[i+1][j] - 2.0* un[i][j] + un[i-1][j]) + (dt/(dy*dy))*(un[i][j+1] - 2.0*un[i][j] + un[i][j-1]));

                v[i][j] = vn[i][j] - un[i][j]*(dt/dx)*(vn[i][j] - vn[i-1][j]) - vn[i][j]*(dt/dy)*(vn[i][j] - vn[i][j-1])
                            - (dt/(rho*2.0*dy)) * (pn[i][j+1] - pn[i][j-1])
                            + nu*((dt/(dx*dx))*(vn[i+1][j] - 2.0* vn[i][j] + vn[i-1][j]) + (dt/(dy*dy))*(vn[i][j+1] - 2.0* vn[i][j] + vn[i][j-1]));
            }
        }

        for (int i = 0; i < nx; i++) {
            u[i][0] = 0.0;
            v[i][0] = 0.0;
            u[i][ny - 1] = 1.0;
            v[i][ny - 1] = 0.0;
        }

        for (int j = 0; j < ny; j++) {
            u[0][j] = 0.0;
            v[0][j] = 0.0;
            u[nx - 1][j] = 0.0;
            v[nx - 1][j] = 0.0;
        }
        
        plt::figure(1); 
        std::vector<double> X1, Y1, P1, U1, V1;

        for (int i = 0; i < ny; ++i) {
            for (int j = 0; j < nx; ++j) {
                X1.push_back(j * dx);
                Y1.push_back(i * dy);
                P1.push_back(p[i][j]);
                U1.push_back(u[j][i]);
                V1.push_back(v[j][i]);
            }
        }

        const long animatedFig = plt::figure(1);
        plt::ion(); 
        plt::clf();
        plt::show();
        plt::quiver(X1, Y1, U1, V1);
        plt::pause(0.1);

    }
    return 0;
}