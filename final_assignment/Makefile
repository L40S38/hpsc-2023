openmp:10_cavity_openmp.cpp
	module purge
	module load gcc
	module load gnuplot
	g++ 10_cavity_openmp.cpp -fopenmp
	./a.out
simd:10_cavity_simd.cpp
	module purge
	module load intel gcc
	module load gnuplot
	g++ 10_cavity_simd.cpp -fopt-info-vec-optimized -march=native -O3
	./a.out
cuda:10_cavity_cuda.cu
	module purge
	module load cuda
	module load gnuplot
	nvcc 10_cavity_cuda.cu -std=c++11
	./a.out
openacc:10_cavity_openacc.cpp
	module purge
	module load nvhpc/22.2
	nvc++ 10_cavity_openacc.cpp -acc
	nvc++ 10_cavity_openacc.cpp -acc -Minfo
	./a.out