//#pragma execution_character_set("utf-8")

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <chrono>

// グローバル変数
const int nx = 41;
const int ny = 41;
const int nt = 500;
//const int nt = 10;
const int nit = 50;
const double dx = 2.0 / (nx - 1);
const double dy = 2.0 / (ny - 1);
const double dt = 0.01;
const double rho = 1.0;
const double nu = 0.02;

double *u, *v, *p, *b, *pn, *un, *vn;

#define u(y,x) u[(y)*ny + (x)]
#define v(y,x) v[(y)*ny + (x)]
#define p(y,x) p[(y)*ny + (x)]
#define b(y,x) b[(y)*ny + (x)]
#define pn(y,x) pn[(y)*ny + (x)]
#define un(y,x) un[(y)*ny + (x)]
#define vn(y,x) vn[(y)*ny + (x)]

#define a(y,x) a[y*ny+x]
#define an(y,x) an[y*ny+x]

// プロットデータをファイルに保存する関数
void sendDataPressure(const std::vector<double>& x, const std::vector<double>& y,
                        //const std::vector<std::vector<double>>& p, 
                        double *p, FILE* fp) {
    //FILE* fp = fopen(filename.c_str(), "w", "utf-8");
    if (fp != nullptr) {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                fprintf(fp, "%f %f %f\n", x[j], y[i], p(i,j));
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
                        //const std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& v,
                        double *u, double *v, FILE* fp) {
    //FILE* fp = fopen(filename.c_str(), "w");
    if (fp != nullptr) {
        for (int i = 1; i < ny; i = i+2) {
            for (int j = 1; j < nx; j = j+2) {
                fprintf(fp, "%f %f %f %f\n", x[j], y[i], u(i,j)*0.5, v(i,j)*0.5);
            }
        }
        fprintf(fp, "e\n");
        //fclose(fp);
    } else {
        std::cerr << "Failed to open file" << std::endl;
    }
}

//境界条件
void border(double *u, double *v) {
    for (int i = 0; i < nx; i++) {
        u(0,i) = 0.0;
        u(ny - 1,i) = 1.0;
        v(0,i) = 0.0;
        v(ny - 1,i) = 0.0;
    }
    for (int i = 0; i < ny; i++) {
        u(i,0) = 0.0;
        u(i,nx - 1) = 0.0;
        v(i,0) = 0.0;
        v(i,nx - 1) = 0.0;
    }
}

// 初期条件の設定
/*
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
}
*/

__global__ void init_zero(double *u, double *v, double *p, double *b){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    u[i] = 0.0;
    v[i] = 0.0;
    p[i] = 0.0;
    b[i] = 0.0;
}

void free_memories(double *u, double *v, double *p, double *b, 
                            double *pn, double *un, double *vn){
    cudaFree(u);
    cudaFree(v);
    cudaFree(p);
    cudaFree(b);
    cudaFree(pn);
    cudaFree(un);
    cudaFree(vn);
}

__global__ void data_insert_b(double *b, double *u, double *v){
    int i = blockIdx.x+1;
    int j = threadIdx.x+1;
    //if(i>0 && i<nx-1 && j>0 && j<ny-1){
    b(j,i) = rho * (1 / dt *
            ((u(j,i+1) - u(j,i-1)) / (2 * dx) + (v(j+1,i) - v(j-1,i)) / (2 * dy)) -
            ((u(j,i+1) - u(j,i-1)) / (2 * dx)) * ((u(j,i+1) - u(j,i-1)) / (2 * dx)) -
            2 * ((u(j+1,i) - u(j-1,i)) / (2 * dy)) *
                ((v(j,i+1) - v(j,i-1)) / (2 * dx)) -
            ((v(j+1,i) - v(j-1,i)) / (2 * dy)) * ((v(j+1,i) - v(j-1,i)) / (2 * dy)));
    //}
}

__global__ void copy_array(double *an, double *a){
    int i = blockIdx.x;
    int j = threadIdx.x;

    an(j,i) = a(j,i);
}

__global__ void data_insert_p(double *p, double *b){
    int i = blockIdx.x+1;
    int j = threadIdx.x+1;
    extern __shared__ double pn[];
    __syncthreads();
    pn(j,i) = p(i);
    __syncthreads();
    //if(i>0 && i<nx-1 && j>0 && j<ny-1){
    p(j,i) = (dy * dy * (pn(j,i+1) + pn(j,i-1)) +
                dx * dx * (pn(j+1,i) + pn(j-1,i)) -
                b(j,i) * dx * dx * dy * dy) /
                (2 * (dx * dx + dy * dy));
    //}
}

__global__ void data_insert_u_v(double *u, double *v, double *un, double *vn, double *p){
    int i = blockIdx.x+1;
    int j = threadIdx.x+1;
    //if(i>0 && i<nx-1 && j>0 && j<ny-1){
    u(j,i) = un(j,i) - un(j,i) * dt / dx * (un(j,i) - un(j,i-1)) -
                                    un(j,i) * dt / dy * (un(j,i) - un(j-1,i)) -
                                    dt / (2 * rho * dx) * (p(j,i+1) - p(j,i-1)) +
                                    nu * dt / (dx * dx) * (un(j,i+1) - 2 * un(j,i) + un(j,i-1)) +
                                    nu * dt / (dy * dy) * (un(j+1,i) - 2 * un(j,i) + un(j-1,i));
    v(j,i) = vn(j,i) - vn(j,i) * dt / dx * (vn(j,i) - vn(j,i-1)) -
                                    vn(j,i) * dt / dy * (vn(j,i) - vn(j-1,i)) -
                                    dt / (2 * rho * dx) * (p(j+1,i) - p(j-1,i)) +
                                    nu * dt / (dx * dx) * (vn(j,i+1) - 2 * vn(j,i) + vn(j,i-1)) +
                                    nu * dt / (dy * dy) * (vn(j+1,i) - 2 * vn(j,i) + vn(j-1,i));
    //}
}
/*
__global__ void data_insert_v(double *v, double *un, double *vn, double *p){
    int i = blockIdx.x+1;
    int j = threadIdx.x+1;
    //if(i>0 && i<nx-1 && j>0 && j<ny-1){
    v(j,i) = vn(j,i) - vn(j,i) * dt / dx * (vn(j,i) - vn(j,i-1)) -
                                    vn(j,i) * dt / dy * (vn(j,i) - vn(j-1,i)) -
                                    dt / (2 * rho * dx) * (p(j+1,i) - p(j-1,i)) +
                                    nu * dt / (dx * dx) * (vn(j,i+1) - 2 * vn(j,i) + vn(j,i-1)) +
                                    nu * dt / (dy * dy) * (vn(j+1,i) - 2 * vn(j,i) + vn(j-1,i));
    //}
}
*/

int main(void){
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
    /*
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx));
    */
    //initialize(u, v, p, b, pn, un, vn);
    cudaMallocManaged(&u, nx*ny*sizeof(double));
    cudaMallocManaged(&v, nx*ny*sizeof(double));
    cudaMallocManaged(&p, nx*ny*sizeof(double));
    cudaMallocManaged(&b, nx*ny*sizeof(double));
    cudaMallocManaged(&pn, nx*ny*sizeof(double));
    cudaMallocManaged(&un, nx*ny*sizeof(double));
    cudaMallocManaged(&vn, nx*ny*sizeof(double));
    init_zero<<<nx,ny>>>(u,v,p,b);
    cudaDeviceSynchronize();
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
        /*
        for(int j=1; j<ny-1; j++){
            for(int i=1; i<nx-1; i++){
                b[j][i] = rho * (1 / dt *
                                 ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
                                 ((u[j][i+1] - u[j][i-1]) / (2 * dx)) * ((u[j][i+1] - u[j][i-1]) / (2 * dx)) -
                                 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy)) *
                                     ((v[j][i+1] - v[j][i-1]) / (2 * dx)) -
                                 ((v[j+1][i] - v[j-1][i]) / (2 * dy)) * ((v[j+1][i] - v[j-1][i]) / (2 * dy)));
            }
        }
        */
        data_insert_b<<<nx-2,ny-2>>>(b,u,v);
        cudaDeviceSynchronize();
        for(int it=0; it<nit; it++){
            //std::vector<std::vector<double>> pn = p;
            //copy_array<<<nx,ny>>>(pn,p);
            //cudaDeviceSynchronize();
            /*
            for(int j=1; j<ny-1; j++){
                for(int i=1; i<nx-1; i++){
                    p[j][i] =
                        (dy * dy * (pn[j][i+1] + pn[j][i-1]) +
                         dx * dx * (pn[j+1][i] + pn[j-1][i]) -
                         b[j][i] * dx * dx * dy * dy) /
                        (2 * (dx * dx + dy * dy));
                }
            }*/
            data_insert_p<<<nx-2,ny-2>>>(p,b);
            cudaDeviceSynchronize();
            for (int i = 0; i < nx; i++) {
                p(0,i) = p(1,i);
                p(ny - 1,i) = 0.0;
            }
            for (int i = 0; i < ny; i++) {
                p(i,0) = p(i,1);
                p(i,nx - 1) = p(i,nx - 2);
            }
        }
        //std::vector<std::vector<double>> un = u;
        //std::vector<std::vector<double>> vn = v;
        copy_array<<<nx,ny>>>(un,u);
        copy_array<<<nx,ny>>>(vn,v);
        cudaDeviceSynchronize();
        /*
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i-1]) -
                                    un[j][i] * dt / dy * (un[j][i] - un[j-1][i]) -
                                    dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1]) +
                                    nu * dt / (dx * dx) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
                                    nu * dt / (dy * dy) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i-1]) -
                                    vn[j][i] * dt / dy * (vn[j][i] - vn[j-1][i]) -
                                    dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i]) +
                                    nu * dt / (dx * dx) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
                                    nu * dt / (dy * dy) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
            }
        }
        */
        data_insert_u_v<<<nx-2,ny-2>>>(u,v,un,vn,p);
        //data_insert_v<<<nx-2,ny-2>>>(v,un,vn,p);
        cudaDeviceSynchronize();
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

    //メモリを解放する
    free_memories(u, v, p, b, pn, un, vn);
    return 0;
}