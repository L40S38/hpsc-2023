#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range,0); 
  for (int i=0; i<n; i++)
#pragma omp atomic update
    bucket[key[i]]++;
  printf("bucket: ");
  for (int i=0;i<range; i++){
    printf("%d ",bucket[i]);
  }
  printf("\n");
  std::vector<int> offset(range,0);
  std::vector<int> tmp(range,0);
  //for (int i=1; i<range; i++) 
  //  offset[i] = offset[i-1] + bucket[i-1];
#pragma omp parallel for
  for (int i=0; i<range; i++){
    offset[i] = bucket[i];
  }

#pragma omp parallel
  for (int j=1; j<range; j<<=1){//Prefix Sumでbucketの0~iまでの合計をoffset[i]に格納
#pragma omp for
    for (int i=0; i<range; i++){
      tmp[i] = offset[i];
    }
#pragma omp for
    for (int i=j; i<range; i++){
      offset[i] += tmp[i-j];
    }
  }
#pragma omp parallel for
  for(int i=0; i<range; i++){
    offset[i] -= bucket[i];//実際に使う値よりもbucket[i]分大きい値が格納されているのでマイナス
  }
  printf("offset: ");
  for (int i=0;i<range; i++){
    printf("%d ",offset[i]);
  }
  printf("\n");

//#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<range; i++) {
    int j = offset[i];
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
