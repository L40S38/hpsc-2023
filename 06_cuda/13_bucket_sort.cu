#include <cstdio>
#include <cstdlib>
#include <vector>
__global__ void zero_init(int *array){/*ゼロ初期化*/
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  array[i] = 0;
  __syncthreads();
}

__global__ void reduction_add(int *array,int *key){/*加算のためのreduction*/
  int i = threadIdx.x;
  int j = blockIdx.x;
  if(key[i]==j)atomicAdd(&array[j],1);
}

__global__ void scan(int *a, int *b, int range) {/*prefix sumのための関数*/
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  for(int j=1; j<range; j<<=1) {
    b[i] = a[i];
    __syncthreads();
    if(i>=j) a[i] += b[i-j];
    __syncthreads();
  }
  b[i] = a[i];
  __syncthreads();
  if(i==0)a[0]=0;
  else a[i]=b[i-1];
}

__global__ void sort(int *offset, int *key, int N) {/*bucket sort*/
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int begin,end;
  if(i==0)begin=0;
  else begin=offset[i-1];
  if(i==N-1)end=N;
  else end=offset[i];

  for(int j=begin;j<end;j++) {
    key[j]=i;
  }
}

int main() {
  int n = 50;
  int range = 5;
  int *key,*bucket,*offset;
  /*メモリ確保*/
  cudaMallocManaged(&key, n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int));
  cudaMallocManaged(&offset, range*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
 
  zero_init<<<1,range>>>(bucket);
  zero_init<<<1,range>>>(offset);
  cudaDeviceSynchronize();
  reduction_add<<<range,n>>>(bucket,key);
  cudaDeviceSynchronize();
  scan<<<1,range>>>(bucket,offset,range);
  cudaDeviceSynchronize();
  /*
  for (int i=0; i<range; i++) {
    printf("%d ",bucket[i]);
  }
  printf("\n");
  for (int i=0; i<range; i++) {
    printf("%d ",offset[i]);
  }
  printf("\n");
  */
  sort<<<1,range>>>(offset,key,n);
  cudaDeviceSynchronize();
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
