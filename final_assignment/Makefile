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
cuda:10_cavity_cuda.cuda
	module purge
	module load cuda
	module load gnuplot
	nvcc 10_cavity_cuda.cu -std=c++11
	./a.out