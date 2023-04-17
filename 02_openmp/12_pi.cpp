#include <cstdio>
#include <omp.h>

int main() {
  int n = 1000;
  double dx = 1. / n;
  double pi = 0;
  double pi_sub[n];

#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<n; i++) {
    double x = (i + 0.5) * dx;
    //pi += 4.0 / (1.0 + x * x) * dx;
    pi_sub[i] = 4.0 / (1.0 + x * x) * dx;
    //printf("%d %d %f\n",omp_get_thread_num(),i,pi_sub[i]);
  }
//#pragma omp parallel
  for (int i=0; i<n; i++) {
    //double x = (i + 0.5) * dx;
    //pi += 4.0 / (1.0 + x * x) * dx;
    pi += pi_sub[i];
    //printf("%d %f\n",i,pi);
  }
  printf("%17.15f\n",pi);
}
