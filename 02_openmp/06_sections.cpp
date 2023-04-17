#include <cstdio>
#include <omp.h>

int main() {
  int i = 1;//外部で初期化されている
#pragma omp parallel num_threads(2)
#pragma omp sections firstprivate(i) //ただのprivateだと、外部で値が設定されていたとしてもゼロ初期化されてしまう
  {//各々でiを引き継いで+1している
#pragma omp section
    printf("%d\n",++i);
#pragma omp section
    printf("%d\n",++i);
  }
}
