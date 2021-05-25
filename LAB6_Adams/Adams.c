#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>

#define Y0 1/(2*log(2))
#define N 768

//Точное решение
double DiffEqSol(double x, double a) {
	return (a + 1) / ((x + 1) * (log(4) - (1 + a) * (log(1 + a)) + (a + 1) * (log(x + 1))));
}

//Дифференциальное уравнение
double DiffEq(double x, double y) { 
	return -y / (x + 1) - y * y;
}


void grid(double* X, double h, double a, double n)
{
	for (int i = 0; i < n; i++)
		X[i] = a + h * i;
}

//Погрешности, возвращает максиммальныую ошибку
double max_error(double* X, double* Y, double a, int n)
{
	double* err = malloc(n * sizeof(double));

	for (int i = 0; i < n; i++)
	{
		err[i] = fabs(DiffEqSol(X[i], a) - Y[i]);
	}

	double max_err = err[0];
	for (int i = 1; i < n; i++)
	{
		if (fabs(err[i]) > max_err)
			max_err = err[i];
	}

	return max_err;
}

void RK3(double* X, double h, int n, double* Y, double a) {

	X[0] = a;
	Y[0] = Y0;

	double K1, K2, K3, K4;

	for (int i = 1; i <= n-1; i++) {
		K1 = h * DiffEq(X[i - 1], Y[i - 1]);
		K2 = h * DiffEq(X[i - 1] + h / 2.0, Y[i - 1] + (h * K1) / 2.0);
		K3 = h * DiffEq(X[i - 1] + h, Y[i - 1] + 2 * h * K2 - h * K1);

		Y[i] = Y[i - 1] + (K1 + 4 * K2 + K3) / 6.0;
	}
}

void Adams3(double* X, double h, int n, double* Y, double a) {

	X[0] = a;
	Y[0] = Y0;

	double* Y_AdamsBashfort = malloc(n * sizeof(double));

	double K1, K2, K3, K4;

	//Вычисляем недостающие для метода Адамса значения y, тк метод трехшаговый
	//нулевой у нас уже есть, соответственно вычислим остальные с таким циклом методом Рунге-Кутты
	for (int i = 1; i <= 2; i++) {
		K1 = h * DiffEq(X[i-1], Y[i-1]);
		K2 = h * DiffEq(X[i-1] + h / 2.0, Y[i-1] + (h * K1) / 2.0);
		K3 = h * DiffEq(X[i-1] + h, Y[i-1] + 2 * h * K2 - h * K1);

		Y[i] = Y[i-1] + (K1 + 4 * K2 + K3) / 6.0;


	}

	//Непосредственно сам Адамс предиктор-корректор
	
	for (int i = 3; i < n; i++)
	{
		Y_AdamsBashfort[i] = Y[i - 1] + (h / 12.0) * (23 * DiffEq(X[i - 1], Y[i - 1]) - 16 * DiffEq(X[i - 2], Y[i - 2]) +
			5 * DiffEq(X[i - 3], Y[i - 3]));
	
		Y[i] = Y[i - 1] + (h / 12.0) * (5 * DiffEq(X[i], Y_AdamsBashfort[i]) + 8 * DiffEq(X[i - 1], Y[i - 1]) - DiffEq(X[i - 2], Y[i - 2]));
	} 

	free(Y_AdamsBashfort);
}

void Adams3RungeRule(double a, double b, double eps) {
	double K1, K2, K3, K4;
	double h = 1; int n = 6;

	double theor_err;
	double fact_err;
	do {
		double* Y_AdamsBashfort = malloc(n * sizeof(double));
		double* TheoretError = malloc(n * sizeof(double));

		double* X = malloc(n * sizeof(double));
		double* Y = malloc(n * sizeof(double));
		Y[0] = Y0;
		grid(X, h, a, n);
		for (int i = 1; i <= 2; i++) {
			K1 = h * DiffEq(X[i - 1], Y[i - 1]);
			K2 = h * DiffEq(X[i - 1] + h / 2.0, Y[i - 1] + (h * K1) / 2.0);
			K3 = h * DiffEq(X[i - 1] + h, Y[i - 1] + 2 * h * K2 - h * K1);
			Y[i] = Y[i - 1] + (K1 + 4 * K2 + K3) / 6.0;
		}

		for (int i = 3; i < n; i++)
		{
			Y_AdamsBashfort[i] = Y[i - 1] + (h / 12.0) * (23 * DiffEq(X[i - 1], Y[i - 1]) - 16 * DiffEq(X[i - 2], Y[i - 2]) +
				5 * DiffEq(X[i - 3], Y[i - 3]));

			Y[i] = Y[i - 1] + (h / 12.0) * (5 * DiffEq(X[i], Y_AdamsBashfort[i]) + 8 * DiffEq(X[i - 1], Y[i - 1]) - DiffEq(X[i - 2], Y[i - 2]));

			TheoretError[i] = fabs(Y_AdamsBashfort[i] - Y[i]) / 10.0;

			theor_err = TheoretError[i];
			if (TheoretError[i] > eps) break;

		}

		if (theor_err < eps)
		{
			FILE* axesX = fopen("X", "w");
			FILE* axesY = fopen("Y-Adams", "w");
			FILE* axesYTrue = fopen("Y-True", "w");

			for (int i = 0; i < n; i++)
				fprintf(axesX, "%.15lf\n", X[i]);
			for (int i = 0; i < n; i++)
				fprintf(axesY, "%.15lf\n", Y[i]);
			for (int i = 0; i < n; i++)
				fprintf(axesYTrue, "%.15lf\n", DiffEqSol(X[i], a));

			fclose(axesY);
			fclose(axesX);
			fclose(axesYTrue);
		}

		printf("THEOR. ERROR = %.15lf \n", theor_err);
		fact_err = max_error(X, Y, a, n);

		n *= 2;
		h = (b - a) / (n - 1);
		free(Y_AdamsBashfort);
		free(Y);
		free(TheoretError);
		free(X);
	} while (theor_err > eps);

	printf("==============\n");
	printf("THEOR. ERROR = %.15lf \n", theor_err);
	printf("FACT. ERROR = %.15lf \n", fact_err);
	printf("n = %d \n", n / 2);

}

void Adams4RungeRule(double a, double b, double eps) {
	double K1, K2, K3, K4;
	double h = 1; int n = 6;

	double theor_err;
	double fact_err;
	do {
		double* Y_AdamsBashfort = malloc(n * sizeof(double));
		double* TheoretError = malloc(n * sizeof(double));

		double* X = malloc(n * sizeof(double));
		double* Y = malloc(n * sizeof(double));
		Y[0] = Y0;
		grid(X, h, a, n);
		for (int i = 1; i <= 3; i++) {
			K1 = h * DiffEq(X[i - 1], Y[i - 1]);
			K2 = h * DiffEq(X[i - 1] + h / 2.0, Y[i - 1] + (h * K1) / 2.0);
			K3 = h * DiffEq(X[i - 1] + h, Y[i - 1] + 2 * h * K2 - h * K1);
			Y[i] = Y[i - 1] + (K1 + 4 * K2 + K3) / 6.0;
		}

		for (int i = 4; i < n; i++)
		{
			Y_AdamsBashfort[i] = Y[i - 1] + (h / 24.0) * (55 * DiffEq(X[i - 1], Y[i - 1]) - 59 * DiffEq(X[i - 2], Y[i - 2]) +
				37 * DiffEq(X[i - 3], Y[i - 3]) - 9 * DiffEq(X[i - 4], Y[i - 4]));

			Y[i] = Y[i - 1] + (h / 24.0) * (9 * DiffEq(X[i], Y_AdamsBashfort[i]) + 19 * DiffEq(X[i - 1], Y[i - 1]) -
				5 * DiffEq(X[i - 2], Y[i - 2]) + DiffEq(X[i - 3], Y[i - 3]));

			TheoretError[i] = (19.0 * fabs(Y_AdamsBashfort[i] - Y[i])) / 270.0;

			theor_err = TheoretError[i];
			if (TheoretError[i] > eps) break;

		}

		if (theor_err < eps)
		{
			FILE* axesX = fopen("X", "w");
			FILE* axesY = fopen("Y-Adams", "w");
			FILE* axesYTrue = fopen("Y-True", "w");

			for (int i = 0; i < n; i++)
				fprintf(axesX, "%.15lf\n", X[i]);
			for (int i = 0; i < n; i++)
				fprintf(axesY, "%.15lf\n", Y[i]);
			for (int i = 0; i < n; i++)
				fprintf(axesYTrue, "%.15lf\n", DiffEqSol(X[i], a));

			fclose(axesY);
			fclose(axesX);
			fclose(axesYTrue);
		}

		printf("THEOR. ERROR = %.15lf \n", theor_err);
		fact_err = max_error(X, Y, a, n);

		n *= 2;
		h = (b - a) / (n - 1);
		free(Y_AdamsBashfort);
		free(Y);
		free(TheoretError);
		free(X);
	} while (theor_err > eps);

	printf("==============\n");
	printf("THEOR. ERROR = %.15lf \n", theor_err);
	printf("FACT. ERROR = %.15lf \n", fact_err);
	printf("n = %d \n", n / 2);

}


void standart_test(double a, double b)
{

	double X[N] = { 0 };
	double Y[N] = { 0 };
	double h;

	h = (b - a) / (N - 1);
	grid(X, h, a, N);

	Adams3(X, h, N, Y, a);

	FILE* axesX = fopen("X", "w");
	FILE* axesY = fopen("Y-Adams", "w");

	for (int i = 0; i < N; i++)
		fprintf(axesX, "%.15lf\n", X[i]);
	for (int i = 0; i < N; i++)
		fprintf(axesY, "%.15lf\n", Y[i]);

	printf("FACT ERR = %lf", max_error(X, Y, a, N));
	fclose(axesY);
	fclose(axesX);
}

void test_h_error_graphic(double a, double b)
{

	FILE* axesH = fopen("h", "w");
	FILE* axesErr = fopen("Error(h)", "w");
	
	double h = 1;
	int n = 6;

	do  {
		double* X = malloc(n * sizeof(double));
		double* Y = malloc(n * sizeof(double));

		grid(X, h, a, n);
		Adams3(X, h, n, Y, a);

		fprintf(axesH, "%.15lf \n", h);
		fprintf(axesErr, "%.15lf \n", max_error(X,Y,a, n));

		free(X);
		free(Y);

		n = n + 5;
		h = (b-a) / (n - 1);

	} while (h > 0.001);

	fclose(axesH);
	fclose(axesErr);
}

int main() {

	double a, b, eps;

	a = 0;
	b = 5;
	eps = 0.01;
	//standart_test(a, b);
	//test_h_error_graphic(a, b);
	Adams4RungeRule(a, b, eps);
	return 0;
}