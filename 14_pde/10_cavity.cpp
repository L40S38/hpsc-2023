#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

void plotData(const vector<vector<double>>& x, const vector<vector<double>>& y,
              const vector<vector<double>>& p, const vector<vector<double>>& u,
              const vector<vector<double>>& v) {
    ofstream dataFile("data.txt");
    if (!dataFile) {
        cerr << "Failed to open data file." << endl;
        return;
    }

    int nx = x[0].size();
    int ny = x.size();

    // Write data to file
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            dataFile << x[i][j] << " " << y[i][j] << " " << p[i][j] << " " << u[i][j] << " " << v[i][j] << endl;
        }
        dataFile << endl;
    }
    dataFile.close();

    // Plot with gnuplot
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe != nullptr) {
        fprintf(gnuplotPipe, "set pm3d\n");
        fprintf(gnuplotPipe, "set dgrid3d\n");
        fprintf(gnuplotPipe, "splot 'data.txt' u 1:2:3 with pm3d title 'Pressure'\n");
        fprintf(gnuplotPipe, "set arrow from 0,0 to 2,0 nohead\n");
        fprintf(gnuplotPipe, "set arrow from 2,0 to 2,2 nohead\n");
        fprintf(gnuplotPipe, "set arrow from 2,2 to 0,2 nohead\n");
        fprintf(gnuplotPipe, "set arrow from 0,2 to 0,0 nohead\n");
        fprintf(gnuplotPipe, "set title 'Velocity Field'\n");
        fprintf(gnuplotPipe, "set xrange [0:2]\n");
        fprintf(gnuplotPipe, "set yrange [0:2]\n");
        fprintf(gnuplotPipe, "set xlabel 'X'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' u 1:2:(sqrt($4*$4+$5*$5)) with vectors title 'Velocity'\n");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe, "exit\n");
        pclose(gnuplotPipe);
    } else {
        cerr << "Failed to open gnuplot pipe." << endl;
    }
}

int main() {
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);
    double dt = 0.01;
    double rho = 1.0;
    double nu = 0.02;

    vector<double> x(nx);
    vector<double> y(ny);
    vector<vector<double>> u(ny, vector<double>(nx));
    vector<vector<double>> v(ny, vector<double>(nx));
    vector<vector<double>> p(ny, vector<double>(nx));
    vector<vector<double>> b(ny, vector<double>(nx));
    vector<vector<double>> X(ny, vector<double>(nx));
    vector<vector<double>> Y(ny, vector<double>(nx));

    // Generate meshgrid
    for (int i = 0; i < nx; ++i) {
        x[i] = i * dx;
    }
    for (int j = 0; j < ny; ++j) {
        y[j] = j * dy;
    }
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            X[j][i] = x[i];
            Y[j][i] = y[j];
        }
    }

    // Main loop
    for (int n = 0; n < nt; ++n) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                b[j][i] = rho * (1 / dt * ((u[j][i + 1] - u[j][i - 1]) / (2 * dx) +
                                           (v[j + 1][i] - v[j - 1][i]) / (2 * dy)) -
                                 ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) * ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) -
                                 2 * ((u[j + 1][i] - u[j - 1][i]) / (2 * dy)) *
                                         ((v[j][i + 1] - v[j][i - 1]) / (2 * dx)) -
                                 ((v[j + 1][i] - v[j - 1][i]) / (2 * dy)) * ((v[j + 1][i] - v[j - 1][i]) / (2 * dy)));
            }
        }
        for (int it = 0; it < nit; ++it) {
            vector<vector<double>> pn = p;
            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {
                    p[j][i] =
                        (dy * dy * (pn[j][i + 1] + pn[j][i - 1]) + dx * dx * (pn[j + 1][i] + pn[j - 1][i]) -
                         b[j][i] * dx * dx * dy * dy) /
                        (2 * (dx * dx + dy * dy));
                }
            }
            for (int i = 0; i < nx; ++i) {
                p[0][i] = p[1][i];
                p[ny - 1][i] = 0.0;
            }
            for (int j = 0; j < ny; ++j) {
                p[j][0] = p[j][1];
                p[j][nx - 1] = p[j][nx - 2];
            }
        }
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                u[j][i] = u[j][i] - u[j][i] * dt / dx * (u[j][i] - u[j][i - 1]) -
                          v[j][i] * dt / dy * (u[j][i] - u[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j][i + 1] - p[j][i - 1]) +
                          nu * dt / dx * dx * (u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) +
                          nu * dt / dy * dy * (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]);
                v[j][i] = v[j][i] - u[j][i] * dt / dx * (v[j][i] - v[j][i - 1]) -
                          v[j][i] * dt / dy * (v[j][i] - v[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j + 1][i] - p[j - 1][i]) +
                          nu * dt / dx * dx * (v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) +
                          nu * dt / dy * dy * (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]);
            }
        }
        for (int i = 0; i < nx; ++i) {
            u[0][i] = 0.0;
            u[ny - 1][i] = 1.0;
            v[0][i] = 0.0;
            v[ny - 1][i] = 0.0;
        }
        for (int j = 0; j < ny; ++j) {
            u[j][0] = 0.0;
            u[j][nx - 1] = 0.0;
            v[j][0] = 0.0;
            v[j][nx - 1] = 0.0;
        }

        // Output data for plotting
        plotData(X, Y, p, u, v);
    }

    return 0;
}
