#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include <cuda.h>
#include <math.h>
#include <curand.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        std::stringstream err;
        err << "GPUassert -- " << file << "(" << line << "): " << cudaGetErrorString(code) << std::endl;
        std::cout << err.str(); // Needed on Windows for some reason
        throw std::runtime_error(err.str());
    }
} 

// The GPU kernel to check if the randomly generated (x,y) pairs fall in the unit circle in the first quadrant or not.
// This is assuming that the total number of (x,y) pairs is an exact multiple of 1024 so that all threads in every
// thread block will be doing work and no out of bounds memory errors occur.
__global__ void checkIn(float *x, float *y, int *counts) {
    int tid = threadIdx.x + blockIdx.x*blockDim.x;
    
    __shared__ int result[1024];
    
    float r = sqrtf(x[tid]*x[tid] + y[tid]*y[tid]);
    result[threadIdx.x] = int(2 - r);
    __syncthreads();
    
    if (threadIdx.x < 31) {
        for (int i = 0; i < 32; ++i) {
            result[threadIdx.x] += result[(threadIdx.x+1)*32 + i];
        }
    }
    __syncthreads();
    
    if (threadIdx.x == 0) {
        for (int i = 1; i < 32; ++i) {
            result[threadIdx.x] += result[i];
        }
        counts[blockIdx.x] = result[threadIdx.x];
    }
}

// This GPU kernel takes the results from checkIn, does a reduction of the results from each block and calculates a
// Monte Carlo estimate of pi which is stored in the array pis.
__global__ void getPi(int *counts, double *pis, int iter) {
    for (int i = 0; i < 32; ++i) {
        counts[threadIdx.x] += counts[(threadIdx.x + 1)*32 + i];
    }
    __syncthreads();
    
    if (threadIdx.x == 0) {
        for (int i = 1; i < 32; ++i) {
            counts[0] += counts[i];
        }
        pis[iter] = 4.0*double(counts[0])/1048576.0;
    }
}

int main() {
    size_t N_blocks = 1024;
    size_t N_threads = 1024;
    int iter = 1000;
    size_t N_pairs = N_blocks*N_threads;
    int *d_counts;
    float *d_x, *d_y;
    std::vector<double> pis(iter);
    std::vector<float> x(N_pairs), y(N_pairs);
    std::vector<int> counts(N_blocks);
    double *d_pis;
    curandGenerator_t gen;
    std::random_device seeder;
    
    cudaMalloc((void**)&d_counts, N_blocks*sizeof(int));
    cudaMalloc((void**)&d_x, N_pairs*sizeof(float));
    cudaMalloc((void**)&d_y, N_pairs*sizeof(float));
    cudaMalloc((void**)&d_pis, iter*sizeof(double));
    
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);
    curandSetPseudoRandomGeneratorSeed(gen, seeder());
    
    double t_1 = omp_get_wtime();
    for (int i = 0; i < iter; ++i) {
        curandGenerateUniform(gen, d_x, N_pairs);
        curandGenerateUniform(gen, d_y, N_pairs);
        
        checkIn<<<N_blocks,N_threads>>>(d_x, d_y, d_counts);
        
        getPi<<<1,31>>>(d_counts, d_pis, i);
    }
    double t_2 = omp_get_wtime();
    std::cout.precision(15);
    std::cout << "Time on GPU: " << t_2 - t_1 << " s\n";
    
    cudaMemcpy(pis.data(), d_pis, iter*sizeof(double), cudaMemcpyDeviceToHost);
    
    std::ofstream fout("pisCUDA.dat");
    fout.precision(15);
    for (int i = 0; i < iter; ++i) {
        //std::cout << pis[i] << "\n";
        fout << pis[i] << "\n";
    }
    fout.close();
    
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_counts);
    cudaFree(d_pis);
    curandDestroyGenerator(gen);
    
    return 0;
}
