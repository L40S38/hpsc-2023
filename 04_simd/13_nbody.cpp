#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], idx[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    idx[i] = i;
  }
  for(int i=0; i<N; i++) {

    //ベクトルの初期化
    __m256 xvec = _mm256_load_ps(x);
    __m256 yvec = _mm256_load_ps(y);
    __m256 mvec = _mm256_load_ps(m);
    __m256 fxvec = _mm256_load_ps(fx);
    __m256 fyvec = _mm256_load_ps(fy);
    __m256 idxvec = _mm256_load_ps(idx);
    __m256 ivec = _mm256_set1_ps(i);
    __m256 zerovec = _mm256_setzero_ps();

    //条件分岐のmask
    __m256 mask = _mm256_cmp_ps(idxvec,ivec,_CMP_NEQ_OQ);

    //rの計算
    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 yivec = _mm256_set1_ps(y[i]);
    __m256 rxvec = _mm256_sub_ps(xivec,xvec);
    __m256 ryvec = _mm256_sub_ps(yivec,yvec);
    __m256 rvec = _mm256_rsqrt_ps(_mm256_add_ps(_mm256_mul_ps(rxvec,rxvec),_mm256_mul_ps(ryvec,ryvec)));
    //rvecはrの逆数
    rvec = _mm256_blendv_ps(zerovec,rvec,mask);

    //代入
    __m256 r3 = _mm256_mul_ps(_mm256_mul_ps(rvec,rvec),rvec);
    __m256 fx_sub = _mm256_mul_ps(_mm256_mul_ps(rxvec,mvec),r3);
    __m256 fy_sub = _mm256_mul_ps(_mm256_mul_ps(ryvec,mvec),r3);
    fx_sub = _mm256_blendv_ps(zerovec,fx_sub,mask);
    fy_sub = _mm256_blendv_ps(zerovec,fy_sub,mask);
    for (int j=1 ; j<N; j<<=1) {//合計値を計算
      if (j==1){
        fx_sub = _mm256_add_ps(fx_sub,_mm256_permute2f128_ps(fx_sub,fx_sub,1));
        fy_sub = _mm256_add_ps(fy_sub,_mm256_permute2f128_ps(fy_sub,fy_sub,1));
      }else{
        fx_sub = _mm256_hadd_ps(fx_sub,fx_sub);
        fy_sub = _mm256_hadd_ps(fy_sub,fy_sub);
      }
    }
    fx_sub = _mm256_blendv_ps(fx_sub,zerovec,mask);
    fy_sub = _mm256_blendv_ps(fy_sub,zerovec,mask);
    fxvec = _mm256_sub_ps(fxvec,fx_sub);
    fyvec = _mm256_sub_ps(fyvec,fy_sub);
    _mm256_store_ps(fx,fxvec);
    _mm256_store_ps(fy,fyvec);
    
    /*
    for(int j=0; j<N; j++) {
      if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
    */
    /*
    0 7.84761 9.13108
    1 60.7831 -40.008
    2 16.7189 -15.641
    3 -20.2919 9.88576
    4 -5.82067 -0.934265
    5 6.15505 -11.1779
    6 6.45807 8.09285
    7 -0.684858 -5.94347
    */
    
    
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
