#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <assert.h>
#include <device_functions.h>
#include <cuda.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define BLOCK_SIZE_X 4
#define BLOCK_SIZE_Y 4
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

__device__ __host__  float q(float x_, float h_x, float y_, float h_y,int counter_for_cells) {
	if (powf(x_, 2) + powf(y_, 2) <= powf(r, 2)) {
		return 1 / (counter_for_cells * h_x * h_y);
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
	
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	if (i > 0 && i < Nx - 1)
	{
		if (j > 0 && j < Ny - 1)
		{
			y_x_layer[j * Nx + i] = (c_x * (t_prev_layer[j * Nx + i + 1] - 2 * t_prev_layer[j * Nx + i] + t_prev_layer[j * Nx + i - 1]) +
				c_y * (t_prev_layer[(j + 1) * Nx + i] - 2 * t_prev_layer[j * Nx + i] + t_prev_layer[(j - 1) * Nx + i]) +
				t_prev_layer[j * Nx + i] + q(x[i], h_x, y[j], h_y, counter_for_cells));
		}
	}

	if (i == 0) {
		y_x_layer[j * Nx + 0] = (gamma_x(y[j], t_i)[0] * h_x - alpha_x(y[j], t_i)[0] * y_x_layer[j * Nx + 1]) / (beta_x(y[j], t_i)[0] * h_x - alpha_x(y[j], t_i)[0]);
	}
	if (i == Nx - 1) {
		y_x_layer[j * Nx + Nx - 1] = (gamma_x(y[j], t_i)[1] * h_x - alpha_x(y[j], t_i)[1] * y_x_layer[j * Nx + Nx - 2]) / (beta_x(y[j], t_i)[1] * h_x - alpha_x(y[j], t_i)[1]);
	}


	if (j == 0) {
		y_x_layer[0 * Nx + i] = (gamma_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0] * y_x_layer[1 * Nx + i]) / \
			(beta_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0]);
	}
	if (j == Ny - 1) {
		y_x_layer[Nx * (Ny - 1) + i] = (gamma_y(x[i], t_i)[1] * h_y - alpha_y(x[i], t_i)[1] * y_x_layer[Nx*(Ny - 2) + i]) / \
			(beta_y(x[i], t_i)[1] * h_y - alpha_y(x[i], t_i)[1]);
	}
}



int main(int argc, char **argv) {
	

	
	float xmin = -1;
	float xmax = 1;
	float ymin = -2;
	float ymax = 2;

	int Nx = 21;
	int Ny = 21;
	float h_x = (xmax-xmin)/(Nx-1);
	float h_y = (ymax-ymin)/(Ny-1);

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
	
	for (int t_i = 0; t_i < Nt; t_i++) {
		solution << <numBlocks, threadsPerBlock >>> (dev_x, Nx, h_x, dev_y, Ny, h_y, dev_t_prev_layer, dev_y_x_layer, t_i, c_x, c_y, counter_for_cells);
		swap(dev_y_x_layer, dev_t_prev_layer);
	}

	
	cudaMemcpy(t_prev_layer, dev_t_prev_layer, mem_size, cudaMemcpyDeviceToHost);
	printf("Reshenie \n");

	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			if (t_prev_layer[j * Nx + i] != 0)
				printf(" j = %d i = %d  t_prev = %f \n", j, i, t_prev_layer[j * Nx + i]);
		}
		printf("\n");
	}
}