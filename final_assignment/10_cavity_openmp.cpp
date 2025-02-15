//#pragma execution_character_set("utf-8")

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
// 時間計測
#include <chrono>
#include <omp.h>

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
// 先に計算しておいた方がいい値リスト
const double rhodx2 = 2 * rho * dx;
const double nudt = nu * dt;
const double dydy = dy * dy;
const double dxdx = dx * dx;


// プロットデータをファイルに保存する関数
void sendDataPressure(const std::vector<double>& x, const std::vector<double>& y,
                        const std::vector<std::vector<double>>& p, FILE* fp) {
    //FILE* fp = fopen(filename.c_str(), "w", "utf-8");
    if (fp != nullptr) {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                fprintf(fp, "%f %f %f\n", x[j], y[i], p[i][j]);
            }
            //fprintf(fp,"\n");
        }
        fprintf(fp, "e\n");
        //fclose(fp);
    } else {
        //std::cerr << "Failed to open file: " << filename << std::endl;
        std::cerr << "Failed to open file" << std::endl;
    }
}

void sendDataVerocity(const std::vector<double>& x, const std::vector<double>& y,
                        const std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& v,
                        FILE* fp) {
    //FILE* fp = fopen(filename.c_str(), "w");
    if (fp != nullptr) {
        for (int i = 1; i < ny; i = i+2) {
            for (int j = 1; j < nx; j = j+2) {
                fprintf(fp, "%f %f %f %f\n", x[j], y[i], u[i][j]*0.5, v[i][j]*0.5);
            }
        }
        fprintf(fp, "e\n");
        //fclose(fp);
    } else {
        std::cerr << "Failed to open file" << std::endl;
    }
}

//境界条件
void border(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v) {
    for (int i = 0; i < nx; i++) {
        u[0][i] = 0.0;
        u[ny - 1][i] = 1.0;
        v[0][i] = 0.0;
        v[ny - 1][i] = 0.0;
    }
    for (int i = 0; i < ny; i++) {
        u[i][0] = 0.0;
        u[i][nx - 1] = 0.0;
        v[i][0] = 0.0;
        v[i][nx - 1] = 0.0;
    }
}

// 初期条件の設定
void initialize(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v,
                std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& b) {
#pragma omp parallel
    for (int i = 0; i < ny; i++) {
#pragma omp for
        for (int j = 0; j < nx; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            b[i][j] = 0.0;
        }
    }
}

int main(void){
    //threadの数の指定
    omp_set_num_threads(10);
    // x軸とy軸
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
    std::chrono::steady_clock::time_point tic, toc;
    double time;

    // gnuplotのパイプラインの作成
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe == nullptr) {
        std::cerr << "Failed to open gnuplot pipe." << std::endl;
        return 1;
    } else {
        // git animate settings
        fprintf(gnuplotPipe, "reset\n");
        //fprintf(gnuplotPipe, "set terminal gif animate\n");
        //fprintf(gnuplotPipe, "set output '10_cavity.gif'\n");
    }

    for(int n=0; n<nt; n++){
        //タイム計測
        toc = std::chrono::steady_clock::now();
#pragma omp parallel
        for(int j=1; j<ny-1; j++){
#pragma omp for schedule(dynamic)
            for(int i=1; i<nx-1; i++){
                b[j][i] = rho * (1 / dt *
                                 ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
                                 ((u[j][i+1] - u[j][i-1]) / (2 * dx)) * ((u[j][i+1] - u[j][i-1]) / (2 * dx)) -
                                 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy)) *
                                     ((v[j][i+1] - v[j][i-1]) / (2 * dx)) -
                                 ((v[j+1][i] - v[j-1][i]) / (2 * dy)) * ((v[j+1][i] - v[j-1][i]) / (2 * dy)));
            }
        }
#pragma omp parallel
        for(int it=0; it<nit; it++){
            std::vector<std::vector<double>> pn = p;
#pragma omp for schedule(dynamic)
            for(int j=1; j<ny-1; j++){
                for(int i=1; i<nx-1; i++){
                    p[j][i] =
                        (dydy * (pn[j][i+1] + pn[j][i-1]) +
                         dxdx * (pn[j+1][i] + pn[j-1][i]) -
                         b[j][i] * dxdx * dydy) /
                        (2 * (dxdx + dydy));
                }
            }
#pragma omp for
            for (int i = 0; i < nx; i++) {
                p[0][i] = p[1][i];
                p[ny - 1][i] = 0.0;
            }
#pragma omp for
            for (int i = 0; i < ny; i++) {
                p[i][0] = p[i][1];
                p[i][nx - 1] = p[i][nx - 2];
            }
        }
        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;
#pragma omp parallel
        for (int j=1; j<ny-1; j++) {
#pragma omp for schedule(dynamic)
            for (int i=1; i<nx-1; i++) {
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i-1]) -
                                    un[j][i] * dt / dy * (un[j][i] - un[j-1][i]) -
                                    dt / (rhodx2) * (p[j][i+1] - p[j][i-1]) +
                                    nudt / (dxdx) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
                                    nudt / (dydy) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i-1]) -
                                    vn[j][i] * dt / dy * (vn[j][i] - vn[j-1][i]) -
                                    dt / (rhodx2) * (p[j+1][i] - p[j-1][i]) +
                                    nudt / (dxdx) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
                                    nudt / (dydy) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
                /*
                2 * rho * dx, nu * dt, dy * dy, dx * dx は先に計算しておいた方が良さそう
                */
            }
        }
        border(u,v);
        // 時間計測, gnuplotの部分は除外
        tic = std::chrono::steady_clock::now();
        time = std::chrono::duration<double>(tic-toc).count();
        std::cout << n + 1 << "," << time << std::endl;

        // gnuplotにプロットコマンドを送信
        fprintf(gnuplotPipe, "set title 'Pressure'\n");
        fprintf(gnuplotPipe, "set xrange [0:2]\n");
        fprintf(gnuplotPipe, "set xlabel 'x'\n");
        fprintf(gnuplotPipe, "set yrange [0:2]\n");
        fprintf(gnuplotPipe, "set ylabel 'y'\n");
        fprintf(gnuplotPipe, "set cbrange [-0.6:0.6]\n");
        fprintf(gnuplotPipe, "set palette maxcolors 24\n");
        fprintf(gnuplotPipe, "set palette defined (-0.6 'blue', 0 'white', 0.6 'red')\n");
        fprintf(gnuplotPipe, "set nokey\n");

        fprintf(gnuplotPipe, "plot '-' u 1:2:3 with image, '-' with vector lc black\n");
        sendDataPressure(x, y, p, gnuplotPipe);
        sendDataVerocity(x, y, u, v, gnuplotPipe);
        fprintf(gnuplotPipe, "pause .01\n");
        fprintf(gnuplotPipe, "\n\n");
        fflush(gnuplotPipe);
    }

    // gnuplotパイプラインを閉じる
    pclose(gnuplotPipe);

    return 0;
}