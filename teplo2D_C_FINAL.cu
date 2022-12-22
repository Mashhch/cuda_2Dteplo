
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;
// условие устойчивости Куранта
float Kurant_condition(float h_x, float h_y, float a) {
	float t_x = powf(h_x, 2) / 2 / powf(a, 2) / 2;
	float t_y = powf(h_y, 2) / 2 / powf(a, 2) / 2;
	float t = (t_x > t_y) ? t_y : t_x;
	return t;
}


// Количество точек, удовлетворяющих условию для круга
float count_for_cells(float* x, int Nx, float* y, int Ny, float r) {
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

float q(float x_, float h_x, float y_, float h_y, float r, int counter_for_cells) {
	if (powf(x_, 2) + powf(y_, 2) <= powf(r, 2)) {
		return 1 / (counter_for_cells * h_x * h_y);
	}
	else
		return 0;
}

float* alpha_x(float x, float y) {
	float alpha[2] = { 1, 1 };
	return alpha;
}

float* beta_x(float x, float y) {
	float beta[2] = { 0, 0 };
	return beta;
}

float* alpha_y(float x, float y) {
	float alpha[2] = { 1, 1 };
	return alpha;
}

float* beta_y(float x, float y) {
	float beta[2] = { 0, 0 };
	return beta;
}

float* gamma_x(float x, float y) {
	float gamma[2] = { 0, 0 };
	return gamma;
}

float* gamma_y(float x, float y) {
	float gamma[2] = { 0, 0 };
	return gamma;
}

float* solution(float* x, int Nx, float h_x, float* y, int Ny, float h_y, float* t_prev_layer, float* y_x_layer, int t_i, float a, float r) {
	float t = Kurant_condition(h_x, h_y, 1);
	int counter_for_cells = count_for_cells(x, Nx, y, Ny, r);
	
	float c_x = powf(a, 2) * t / powf(h_x, 2);
	float c_y = powf(a, 2) * t / powf(h_y, 2);


	for (int i = 1; i < Nx - 1; i++) {
		for (int j = 1; j < Ny - 1; j++) {
			y_x_layer[j * Nx + i] = (c_x * (t_prev_layer[j * Nx + i + 1] - 2 * t_prev_layer[j * Nx + i] + t_prev_layer[j * Nx + i - 1]) +
				c_y * (t_prev_layer[(j + 1) * Nx + i] - 2 * t_prev_layer[j * Nx + i] + t_prev_layer[(j - 1) * Nx + i]) +
				t_prev_layer[j * Nx + i] + q(x[i], h_x, y[j], h_y, r, counter_for_cells));
		}
	}


	for (int j = 1; j < Ny - 1; j++) {
		y_x_layer[j * Nx + 0] = (gamma_x(y[j], t_i)[0] * h_x - alpha_x(y[j], t_i)[0] * y_x_layer[j * Nx + 1]) /(beta_x(y[j], t_i)[0] * h_x - alpha_x(y[j], t_i)[0]);
		y_x_layer[j * Nx + Nx - 1] = (gamma_x(y[j], t_i)[1] * h_x - alpha_x(y[j], t_i)[1] * y_x_layer[j * Nx + Nx - 2]) /(beta_x(y[j], t_i)[1] * h_x - alpha_x(y[j], t_i)[1]);		
	}


	for (int i = 0; i < Nx; i++) {
		y_x_layer[0 * Nx + i] = (gamma_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0] * y_x_layer[1 * Nx + i]) /\
			(beta_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0]);
		y_x_layer[Nx * (Ny - 1) + i] = (gamma_y(x[i], t_i)[1] * h_y - alpha_y(x[i], t_i)[1] * y_x_layer[Nx*(Ny - 2) + i]) / \
			(beta_y(x[i], t_i)[1] * h_y - alpha_y(x[i], t_i)[1]);
	}


	
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			if (y_x_layer[j * Nx + i] != 0)
				printf("time t = %d j = %d i = %d  u = %f \n", t_i, j, i, y_x_layer[j * Nx + i]);
		}
	}
	return y_x_layer;
}



int main(int argc, char **argv) {

	float* my;
	float a = 1;
	float r = 0.25;
	float xmin = -1;
	float xmax = 1;
	float ymin = -2;
	float ymax = 2;
	float h_x = 0.1;
	float h_y = 0.2;
	int Nt = 10;
	int Nx = 21;
	int Ny = 21;
	float* x = (float*)malloc(sizeof(float) * Nx);
	float* y = (float*)malloc(sizeof(float) * Ny);
	float* t = (float*)malloc(sizeof(float) * Nt);
	for (int i = 0; i < Nx; i++) {
		x[i] = xmin + h_x * i;
	}
	for (int j = 0; j < Ny; j++) {
		y[j] = ymin + h_y * j;
	}
	float* y_x_layer = (float*)malloc(Nx * Ny * sizeof(float));
	float* t_prev_layer = (float*)malloc(Nx * Ny * sizeof(float));
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			y_x_layer[j * Nx + i] = 0;
			t_prev_layer[j * Nx + i] = 0;
		}
	}
	float* swap;
	for (int t_i = 0; t_i < Nt; t_i++) {
		swap = solution(x, Nx, h_x, y, Ny, h_y, t_prev_layer, y_x_layer, t_i, a, r);
		y_x_layer = t_prev_layer;
		t_prev_layer = swap;

	}

	printf("Reshenie \n");
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			if (t_prev_layer[j * Nx + i] != 0)
				printf(" j = %d i = %d  t_prev = %f \n",  j, i, t_prev_layer[j * Nx + i]);
		}
		printf("\n");
	}
}