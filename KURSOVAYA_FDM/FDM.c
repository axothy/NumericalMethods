#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define FALSE 0
#define TRUE 1
double p(double x) {
	return 1;
}

double q(double x) {
	return (-4 * (x * x + 3.00)) / (x * x * x + 6 * x);
}

double r(double x) {
	return (7 * x) / (x * x * x + 6 * x);
}

double f(double x) {
	return (x * x * x * x) / (x * x * x + 6 * x);
}

double DiffEqSol(double x) {
	return x * x * x;
}

double GetMaxError(double* x, double* y, int n) {
	double errNorm = 0;
	double* error = malloc(n * sizeof(double));

	for (int i = 0; i < n; i++) {
		error[i] = (fabs(y[i] - DiffEqSol(x[i]))) / (DiffEqSol(x[i]));
	}
	errNorm = error[0];
	for (int i = 1 ; i < n ; i++)
	{
		if (error[i] > errNorm)
			errNorm = error[i];
	}

	FILE* Errors = fopen("Errors", "w");

	for (int i = 0; i < n; i++)
	{
		fprintf(Errors, "%.15lf\n", error[i]);
	}
	free(error);
	return errNorm;
}


int FiniteDifferenceMethod(double (*p)(double x), double (*q)(double x), double (*r)(double x), double (*f)(double x), double a, double b, double y0, double yn, double h, double** y, double** x, int n) {
	int k = n;
	double* F, * C, * D, * B, * alpha, * beta;
	double pk, qk;

	//выделяем память под векторы и матрицу системы
	*x = (double*)malloc(sizeof(double) * (n + 1));
	*y = (double*)malloc(sizeof(double) * (n + 1));
	F = (double*)malloc(sizeof(double) * (n + 1));
	C = (double*)malloc(sizeof(double) * (n + 1));
	D = (double*)malloc(sizeof(double) * (n + 1));
	B = (double*)malloc(sizeof(double) * (n + 1));
	alpha = (double*)malloc(sizeof(double) * (n + 1));
	beta = (double*)malloc(sizeof(double) * (n + 1));


	//заполняем трехдиагональную матрицу коэффициентов системы
	B[0] = 0;
	C[0] = 1;
	D[0] = 0;
	F[0] = y0;
	(*x)[0] = a;
	for (k = 1; k < n; k++) {
		(*x)[k] = (*x)[k - 1] + h;

		pk = p((*x)[k]);
		qk = q((*x)[k]);
		B[k] = pk - h * qk / 2.0;
		C[k] = r((*x)[k]) * h * h - 2.0 * pk;
		D[k] = pk + h * qk / 2.0;

		F[k] = f((*x)[k]) * h * h;
	}
	B[n] = 0;
	C[n] = 1;
	D[n] = 0;
	F[n] = yn;
	(*x)[n] = b;

	//решаем систему методом прогонки
	//прямой ход
	beta[0] = -D[0] / C[0];
	alpha[0] = F[0] / C[0];
	for (k = 1; k < n; k++) {
		beta[k] = -D[k] / (B[k] * beta[k - 1] + C[k]);
		alpha[k] = (F[k] - B[k] * alpha[k - 1]) / (B[k] * beta[k - 1] + C[k]);
	}
	beta[n] = 0;
	alpha[n] = (F[n] - B[n] * alpha[n - 1]) / (B[n] * beta[n - 1] + C[n]);
	//обратный ход
	(*y)[n] = alpha[n];
	for (k = n - 1; k >= 0; k--)
		(*y)[k] = beta[k] * (*y)[k + 1] + alpha[k];

	free(F);
	free(C);
	free(D);
	free(B);
	free(alpha);
	free(beta);
	return TRUE;
}

int get_theor_err(double eps) {
	FILE* axesY1 = fopen("y1.txt", "r");
	FILE* axesY2 = fopen("y2.txt", "r");

	double y1;
	fscanf(axesY1, "%lf", &y1);

	double y2;
	fscanf(axesY2, "%lf", &y2);

	double err = 0;
	double err1 = 0;

	while (fscanf(axesY2, "%lf", &y2) == 1) {
		fscanf(axesY1, "%lf", &y1);
		fscanf(axesY2, "%lf", &y2);

		err1 = fabs(y1 - y2) / 3;
		if (err1 > eps) {
			fclose(axesY1);
			fclose(axesY2);
			return 1;
		}
	}

	fclose(axesY1);
	fclose(axesY2);

	return 0;
}

int PrintSolutionInFile(double* x, double* y, int n, const char* fileX, const char* fileY) {
	FILE* axesX = fopen(fileX, "w");
	FILE* axesY = fopen(fileY, "w");
	for (int i = 0; i < n; i++)
	{
		fprintf(axesX, "%.15lf\n", x[i]);
		fprintf(axesY, "%.15lf\n", y[i]);
	}
	return TRUE;
}

void test_h_error_graphic(double a, double b, double y0, double yn)
{
	double* x = NULL, * y = NULL;
	FILE* axesH = fopen("h", "w");
	FILE* axesErr = fopen("Error(h)", "w");

	double h=0.2;
	int n = 6;
	do {
		FiniteDifferenceMethod(p, q, r, f, a, b, y0, yn, h, &y, &x, n);

		fprintf(axesH, "%.15lf \n", h);
		fprintf(axesErr, "%.15lf \n", GetMaxError(x, y, n));

		free(x);
		free(y);

		n = n + 5;
		h = (b - a) / (n - 1);
	} while (h > 0.001);

	fclose(axesH);
	fclose(axesErr);
}


int main() {
	double a = 0.1, b = 1.0;
	double y0 = 0.00101, yn = 1.01, h;
	double* x = NULL, * y = NULL;
	int n;

	h = 0.01;
	n = (int)((b - a) / h);
	FiniteDifferenceMethod(p, q, r, f, a, b, y0, yn, h, &y, &x, n);
	PrintSolutionInFile(x, y, n, "X.txt", "Y.txt");
	printf("MAX ERR = %.15lf", GetMaxError(x, y, n));

	//test_h_error_graphic(a, b, y0, yn);

	return 0;
}