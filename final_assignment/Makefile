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