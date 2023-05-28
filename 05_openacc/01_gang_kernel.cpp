#include <cstdio>
#include <openacc.h>

int main() {
#pragma acc kernels
  for(int i=0; i<8; i++) {
    printf("(%d,%d,%d): %d\n",
           __pgi_gangidx(),
           __pgi_workeridx(),
           __pgi_vectoridx(),i);
  }
}
