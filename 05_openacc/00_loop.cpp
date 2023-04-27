#include <cstdio>
#include <openacc.h>

int main() {
#pragma acc parallel loop vector
  for(int i=0; i<8; i++) {
    printf("%d: %d\n",__pgi_vectoridx(),i);
  }
}
