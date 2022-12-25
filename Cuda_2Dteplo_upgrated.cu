#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <assert.h>
#include <device_functions.h>
#include <cuda.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>  

#define BLOCK_SIZE_X 32
#define BLOCK_SIZE_Y 32
#define a 1
#define Nt 10
#define r 0.25



__device__ __host__ float Kurant_condition(float h_x, float h_y) {
	float t_x = powf(h_x, 2) / 2 / powf(a, 2) / 2;
	float t_y = powf(h_y, 2) / 2 / powf(a, 2) / 2;
	float t = (t_x > t_y) ? t_y : t_x;
	return t;
}



// Количество точек, удовлетворяющих условию для круга
__device__ __host__  float count_for_cells(float* x, int Nx, float* y, int Ny) {
	int counter_for_cells = 0;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (powf(x[i], 2) + powf(y[j], 2) <= powf(r, 2)) {
				counter_for_cells += 1;
			}
		}
	}
	return counter_for_cells;
}

__device__ __host__  float q(float x_, float h_x, float y_, float h_y, int counter_for_cells) {
	if (powf(x_, 2) + powf(y_, 2) <= powf(r, 2)) {
		return 1 / ((float)counter_for_cells * h_x * h_y);
	}
	else
		return 0;
}

__device__ __host__  float* alpha_x(float x, float y) {
	float alpha[2] = { 1, 1 };
	return alpha;
}

__device__ __host__  float* beta_x(float x, float y) {
	float beta[2] = { 0, 0 };
	return beta;
}

__device__ __host__  float* alpha_y(float x, float y) {
	float alpha[2] = { 1, 1 };
	return alpha;
}

__device__ __host__  float* beta_y(float x, float y) {
	float beta[2] = { 0, 0 };
	return beta;
}

__device__ __host__  float* gamma_x(float x, float y) {
	float gamma[2] = { 0, 0 };
	return gamma;
}

__device__ __host__  float* gamma_y(float x, float y) {
	float gamma[2] = { 0, 0 };
	return gamma;
}

__device__ __host__ void swap(float* &c, float* &b) {
	float *temp = c;
	c = b;
	b = temp;
}

__global__ void solution(float* x, int Nx, float h_x, float* y, int Ny, float h_y, float* t_prev_layer, float* y_x_layer, int t_i, float c_x, float c_y, float counter_for_cells) {

	__shared__ float shared_block[BLOCK_SIZE_X + 2][BLOCK_SIZE_Y + 32];
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int shared_i = threadIdx.x + 1;
	int shared_j = threadIdx.y + 1;

	shared_block[shared_j][shared_i] = t_prev_layer[j * Nx + i];

	if (threadIdx.x == 0 && j != 0) {
		shared_block[shared_j][0] = t_prev_layer[j * Nx + i - 1];
	}
	if (threadIdx.x == BLOCK_SIZE_X - 1 && j != Ny - 1) {
		shared_block[shared_j][BLOCK_SIZE_X + 1] = t_prev_layer[j * Nx + i + 1];
	}
	if (threadIdx.y == 0 && i != 0) {
		shared_block[0][shared_i] = t_prev_layer[(j - 1) * Nx + i];
	}
	if (threadIdx.y == BLOCK_SIZE_Y - 1 && i != Nx - 1) {
		shared_block[BLOCK_SIZE_Y + 1][shared_i] = t_prev_layer[(j + 1) * Nx + i];
	}

	__syncthreads();

	if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1) {
		y_x_layer[j * Nx + i] = (c_x * (shared_block[shared_j][shared_i + 1] - 2 * shared_block[shared_j][shared_i] + shared_block[shared_j][shared_i - 1]) +
			c_y * (shared_block[shared_j + 1][shared_i] - 2 * shared_block[shared_j][shared_i] + shared_block[shared_j - 1][shared_i]) +
			shared_block[shared_j][shared_i] + q(x[i], h_x, y[j], h_y, counter_for_cells));
	}

	__syncthreads();

	if (i == 0) {
		y_x_layer[j * Nx + 0] = (gamma_x(y[j], t_i)[0] * h_x - alpha_x(y[j], t_i)[0] * y_x_layer[j * Nx + 1]) / (beta_x(y[j], t_i)[0] * h_x - alpha_x(y[j], t_i)[0]);
	}
	else if (i == Nx - 1) {
		y_x_layer[j * Nx + Nx - 1] = (gamma_x(y[j], t_i)[1] * h_x - alpha_x(y[j], t_i)[1] * y_x_layer[j * Nx + Nx - 2]) / (beta_x(y[j], t_i)[1] * h_x - alpha_x(y[j], t_i)[1]);
	}


	else if (j == 0) {
		y_x_layer[0 * Nx + i] = (gamma_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0] * y_x_layer[1 * Nx + i]) / \
			(beta_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0]);
	}
	else if (j == Ny - 1) {
		y_x_layer[(Ny - 1)*Nx + i] = (gamma_y(x[i], t_i)[1] * h_y - alpha_y(x[i], t_i)[1] * y_x_layer[(Ny - 2)*Nx + i]) / \
			(beta_y(x[i], t_i)[1] * h_y - alpha_y(x[i], t_i)[1]);
	}


}



int main(int argc, char **argv) {



	float xmin = -1;
	float xmax = 1;
	float ymin = -2;
	float ymax = 2;

	int Nx = 401;
	int Ny = 401;
	float h_x = (xmax - xmin) / (Nx - 1);
	float h_y = (ymax - ymin) / (Ny - 1);

	float t = Kurant_condition(h_x, h_y);
	float c_x = powf(a, 2) * t / powf(h_x, 2);
	float c_y = powf(a, 2) * t / powf(h_y, 2);

	float* x = (float*)malloc(sizeof(float) * Nx);
	float* y = (float*)malloc(sizeof(float) * Ny);


	for (int i = 0; i < Nx; i++) {
		x[i] = xmin + h_x * i;
	}

	for (int j = 0; j < Ny; j++) {
		y[j] = ymin + h_y * j;
	}
	int counter_for_cells = count_for_cells(x, Nx, y, Ny);


	unsigned int mem_size = sizeof(float)*Nx * Ny;
	float* y_x_layer = (float*)malloc(mem_size);
	float* t_prev_layer = (float*)malloc(mem_size);
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			y_x_layer[j * Nx + i] = 0;
			t_prev_layer[j * Nx + i] = 0;
		}
	}

	float *dev_y_x_layer, *dev_t_prev_layer, *dev_x, *dev_y;
	cudaMalloc((void**)&dev_x, Nx * sizeof(float));
	cudaMalloc((void**)&dev_y, Ny * sizeof(float));
	cudaMalloc((void**)&dev_y_x_layer, mem_size);
	cudaMalloc((void**)&dev_t_prev_layer, mem_size);


	cudaMemcpy(dev_x, x, Nx * sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_y, y, Ny * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y_x_layer, y_x_layer, mem_size, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t_prev_layer, t_prev_layer, mem_size, cudaMemcpyHostToDevice);

	dim3 numBlocks(Nx / BLOCK_SIZE_X + 1, Ny / BLOCK_SIZE_Y + 1);
	dim3 threadsPerBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);

	//T_step <<<1,1>>> (dev_x, Nx, h_x, dev_y, Ny, h_y, dev_t_prev_layer, dev_y_x_layer, c_x, c_y, counter_for_cells);
	cudaEvent_t start, stop;
	float time;
	double time_spent = 0.0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	clock_t begin1 = clock();
	for (int t_i = 0; t_i < Nt; t_i++) {
		solution << <numBlocks, threadsPerBlock >> > (dev_x, Nx, h_x, dev_y, Ny, h_y, dev_t_prev_layer, dev_y_x_layer, t_i, c_x, c_y, counter_for_cells);
		swap(dev_y_x_layer, dev_t_prev_layer);
	}

	cudaMemcpy(t_prev_layer, dev_t_prev_layer, mem_size, cudaMemcpyDeviceToHost);
	cudaEventRecord(stop, 0);
	clock_t end1 = clock();

	// рассчитать прошедшее время, найдя разницу (end - begin) и
	// деление разницы на CLOCKS_PER_SEC для перевода в секунды
	time_spent += (double)(end1 - begin1) / CLOCKS_PER_SEC;

	printf("The elapsed time is %f seconds \n", time_spent);
	
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	printf("time = %f \n", time);


	free(x);
	free(y);
	free(t_prev_layer);
	free(y_x_layer);
	cudaFree(dev_t_prev_layer);
	cudaFree(dev_y_x_layer);
	cudaFree(dev_x);
	cudaFree(dev_y);
}