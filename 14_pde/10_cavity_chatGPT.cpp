#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

// グローバル変数
const int nx = 41;
const int ny = 41;
const int nt = 500;
const int nit = 50;
const double dx = 2.0 / (nx - 1);
const double dy = 2.0 / (ny - 1);
const double dt = 0.01;
const double rho = 1.0;
const double nu = 0.02;

// 初期条件の設定
void initialize(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v,
                std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& b) {
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            b[i][j] = 0.0;
        }
    }

    // 境界条件の設定
    for (int i = 0; i < nx; i++) {
        u[0][i] = 0.0;
        u[ny - 1][i] = 0.0;
        v[0][i] = 0.0;
        v[ny - 1][i] = 0.0;
    }
    for (int i = 0; i < ny; i++) {
        u[i][0] = 0.0;
        u[i][nx - 1] = 1.0;
        v[i][0] = 0.0;
        v[i][nx - 1] = 0.0;
    }
}

// プロットデータをファイルに保存する関数
void saveData(const std::vector<std::vector<double>>& data, const std::string& filename) {
    FILE* fp = fopen(filename.c_str(), "w");
    if (fp != nullptr) {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                fprintf(fp, "%f ", data[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    } else {
        std::cerr << "Failed to open file: " << filename << std::endl;
    }
}

int main() {
    // グリッドの作成
    std::vector<double> x(nx);
    std::vector<double> y(ny);
    for (int i = 0; i < nx; i++) {
        x[i] = i * dx;
    }
    for (int i = 0; i < ny; i++) {
        y[i] = i * dy;
    }

    // 変数の初期化
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx));
    initialize(u, v, p, b);

    // gnuplotのパイプラインの作成
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe == nullptr) {
        std::cerr << "Failed to open gnuplot pipe." << std::endl;
        return 1;
    }

    // プロットデータをファイルに保存しながらgnuplotに送信
    for (int n = 0; n < nt; n++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                b[j][i] = rho * (1 / dt *
                                 ((u[j][i + 1] - u[j][i - 1]) / (2 * dx) +
                                  (v[j + 1][i] - v[j - 1][i]) / (2 * dy)) -
                                 ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) *
                                     ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) -
                                 2 * ((u[j + 1][i] - u[j - 1][i]) / (2 * dy)) *
                                     ((v[j][i + 1] - v[j][i - 1]) / (2 * dx)) -
                                 ((v[j + 1][i] - v[j - 1][i]) / (2 * dy)) *
                                     ((v[j + 1][i] - v[j - 1][i]) / (2 * dy)));
            }
        }

        for (int it = 0; it < nit; it++) {
            std::vector<std::vector<double>> pn = p;
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    p[j][i] =
                        (dy * dy * (pn[j][i + 1] + pn[j][i - 1]) +
                         dx * dx * (pn[j + 1][i] + pn[j - 1][i]) -
                         b[j][i] * dx * dx * dy * dy) /
                        (2 * (dx * dx + dy * dy));
                }
            }
            // 境界条件の設定
            for (int i = 0; i < nx; i++) {
                p[0][i] = p[1][i];
                p[ny - 1][i] = 0.0;
            }
            for (int i = 0; i < ny; i++) {
                p[i][0] = p[i][1];
                p[i][nx - 1] = p[i][nx - 2];
            }
        }

        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                u[j][i] = un[j][i] -
                          un[j][i] * dt / dx * (un[j][i] - un[j][i - 1]) -
                          vn[j][i] * dt / dy * (un[j][i] - un[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j][i + 1] - p[j][i - 1]) +
                          nu * (dt / (dx * dx) *
                                    (un[j][i + 1] - 2 * un[j][i] + un[j][i - 1]) +
                                dt / (dy * dy) *
                                    (un[j + 1][i] - 2 * un[j][i] + un[j - 1][i]));

                v[j][i] = vn[j][i] -
                          un[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1]) -
                          vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i]) -
                          dt / (2 * rho * dy) * (p[j + 1][i] - p[j - 1][i]) +
                          nu * (dt / (dx * dx) *
                                    (vn[j][i + 1] - 2 * vn[j][i] + vn[j][i - 1]) +
                                dt / (dy * dy) *
                                    (vn[j + 1][i] - 2 * vn[j][i] + vn[j - 1][i]));
            }
        }

        // プロットデータをファイルに保存
        saveData(p, "pressure.dat");
        saveData(u, "velocity_u.dat");
        saveData(v, "velocity_v.dat");

        // gnuplotにプロットコマンドを送信
        fprintf(gnuplotPipe, "set title 'Pressure'\n");
        fprintf(gnuplotPipe, "set xrange [0:2]\n");
        fprintf(gnuplotPipe, "set yrange [0:2]\n");
        fprintf(gnuplotPipe, "plot 'pressure.dat' matrix with image\n");
        fflush(gnuplotPipe);

        // 一時停止
        std::cout << "Step: " << n + 1 << " / " << nt << std::endl;
    }

    // gnuplotパイプラインを閉じる
    pclose(gnuplotPipe);

    return 0;
}
