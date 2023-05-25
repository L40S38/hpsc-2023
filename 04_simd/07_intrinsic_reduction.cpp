#include <cstdio>
#include <immintrin.h>

int main() {
  const int N = 8;
  float a[N];
  for (int i=0; i<N; i++)
    a[i] = 1;
  __m256 avec = _mm256_load_ps(a);
  __m256 bvec = _mm256_permute2f128_ps(avec,avec,1);
  /*
  bvec = _mm256_add_ps(bvec,avec);
  _mm256_store_ps(a, bvec);
  for (int i=0; i<N; i++)
    printf("%g\n",a[i]);
  bvec = _mm256_hadd_ps(bvec,bvec);
  _mm256_store_ps(a, bvec);
  for (int i=0; i<N; i++)
    printf("%g\n",a[i]);
  bvec = _mm256_hadd_ps(bvec,bvec);
  _mm256_store_ps(a, bvec);
  for (int i=0; i<N; i++)
    printf("%g\n",a[i]);
  */
  for (int j=1; j<N; j<<=1)
    bvec = _mm256_hadd_ps(bvec,bvec);
    _mm256_store_ps(a, bvec);
    for (int i=0; i<N; i++)
      printf("%g\n",a[i]);

}
