#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(){
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);
    double dt = 0.01;
    double rho = 1.0;
    double nu = 0.02;

    std::vector<double> x(nx, 0.0);
    std::vector<double> y(ny, 0.0);

    for (int i = 0; i < nx; i++) {
        x[i] = i * dx;
    }

    for (int i = 0; i < ny; i++) {
        y[i] = i * dy;
    }

    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> X(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> Y(ny, std::vector<double>(nx, 0.0));

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            X[i][j] = x[j];
            Y[i][j] = y[i];
        }
    }

    for (int n = 0; n < nt; n++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                b[j][i] = rho * (1 / dt *
                    ((u[j][i + 1] - u[j][i - 1]) / (2 * dx) + (v[j + 1][i] - v[j - 1][i]) / (2 * dy)) -
                    ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) * ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) -
                    2 * ((u[j + 1][i] - u[j - 1][i]) / (2 * dy)) * ((v[j][i + 1] - v[j][i - 1]) / (2 * dx)) -
                    ((v[j + 1][i] - v[j - 1][i]) / (2 * dy)) * ((v[j + 1][i] - v[j - 1][i]) / (2 * dy))
                );
            }
        }

        for (int it = 0; it < nit; it++) {
            std::vector<std::vector<double>> pn = p;
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    p[j][i] = (dy * dy * (pn[j][i + 1] + pn[j][i - 1]) +
                               dx * dx * (pn[j + 1][i] + pn[j - 1][i]) -
                               b[j][i] * dx * dx * dy * dy) /
                              (2 * (dx * dx + dy * dy));
                }
            }

            for (int i = 0; i < ny; i++) {
                p[i][nx - 1] = p[i][nx - 2];
                p[i][0] = p[i][1];
            }

            for (int i = 0; i < nx; i++) {
                p[0][i] = p[1][i];
                p[ny - 1][i] = 0;
            }

            std::vector<std::vector<double>> un = u;
            std::vector<std::vector<double>> vn = v;

            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    u[j][i] = un[j][i] -
                               un[j][i] * dt / dx * (un[j][i] - un[j][i - 1]) -
                               vn[j][i] * dt / dy * (un[j][i] - un[j - 1][i]) -
                               dt / (2 * rho * dx) * (p[j][i + 1] - p[j][i - 1]) +
                               nu * dt / (dx * dx) * (un[j][i + 1] - 2 * un[j][i] + un[j][i - 1]) +
                               nu * dt / (dy * dy) * (un[j + 1][i] - 2 * un[j][i] + un[j - 1][i]);

                    v[j][i] = vn[j][i] -
                               un[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1]) -
                               vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i]) -
                               dt / (2 * rho * dy) * (p[j + 1][i] - p[j - 1][i]) +
                               nu * dt / (dx * dx) * (vn[j][i + 1] - 2 * vn[j][i] + vn[j][i - 1]) +
                               nu * dt / (dy * dy) * (vn[j + 1][i] - 2 * vn[j][i] + vn[j - 1][i]);
                }
            }

            for (int i = 0; i < nx; i++) {
                u[0][i] = 0;
                u[ny - 1][i] = 1;
                v[0][i] = 0;
                v[ny - 1][i] = 0;
            }

            for (int i = 0; i < ny; i++) {
                u[i][0] = 0;
                u[i][nx - 1] = 0;
                v[i][0] = 0;
                v[i][nx - 1] = 0;
            }

            //plt::contourf(X, Y, p, {{"alpha", 0.5}, {"cmap", "coolwarm"}});
            plt::contour(X, Y, p);
            plt::quiver(X[0], Y[0], u[0], v[0]);
            plt::pause(0.01);
            plt::clf();
        }
        plt::show();
    }

    return 0;
}
