#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>

#define Y0 1/(2*log(2))

//Дифференциальное уравнение
double DiffEqFunc(double x, double y) {
	return -y / (x + 1) - y * y;
}

double DiffEqSol(double x) {
	return 1 / ((1 + x) * log(fabs(1 + x)));
}

//Рунге-Кутта
void RK3(uint64_t n, double a, double b, const char* fileX, const char* fileY)
{
	uint64_t i = 0;

	double h = (b - a) / n;

	double x1 = a;
	double y1 = Y0;
	double x2;
	double y2;

	FILE* axesX = fopen(fileX, "w");
	FILE* axesY = fopen(fileY, "w");

	fprintf(axesX, "%.15lf\n", x1);
	fprintf(axesY, "%.15lf\n", y1);

	double K1, K2, K3;
	for (i = 1; x1 < b; i++) {
		K1 = h * DiffEqFunc(x1, y1);
		K2 = h * DiffEqFunc(x1 + h / 2.0, y1 + (h * K1) / 2.0);
		K3 = h * DiffEqFunc(x1 + h, y1 + 2 * h * K2 - h * K1);
		x2 = x1 + h;
		y2 = y1 + (K1 + 4 * K2 + K3) / 6.0;
		fprintf(axesX, "%.15lf\n", x2);
		fprintf(axesY, "%.15lf\n", y2);
		x1 = x2;
		y1 = y2;
	}
	fclose(axesX);
	fclose(axesY);
}


int get_err(double eps) {
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

		err1 = fabs(y1 - y2) / 7;
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

double get_norm() {
	FILE* axesX2 = fopen("x2.txt", "r");
	FILE* axesY2 = fopen("y2.txt", "r");

	double x2;
	double y2;
	double norm = 0, difference;

	while (fscanf(axesY2, "%lf", &y2) == 1) {
		fscanf(axesX2, "%lf", &x2);

		difference = fabs(DiffEqSol(x2) - y2);
		if (difference > norm) {
			norm = difference;
		}
	}

	return norm;
}

int main()
{
	uint64_t n = 2;
	double a = 1.0, b = 5.0;
	double err = 0;
	double eps = 0.1;

	FILE* EPS = fopen("Epsilon.txt", "w");
	FILE* DEGREE2 = fopen("N.txt", "w");
	FILE* FACTERR = fopen("Fact error.txt", "w");
	int i = 1;

	for (; eps > 1e-6; eps /= 10) {
		do {
			RK3(n, a, b, "x1.txt", "y1.txt");
			n *= 2;
			RK3(n, a, b, "x2.txt", "y2.txt");
			n *= 2;
		} while (get_err(eps));


		fprintf(EPS, "%.15lf\n", eps);
		fprintf(DEGREE2, "%i\n", (int)log2((double)n));
		fprintf(FACTERR, "%.15lf\n", get_norm());

		printf("%i. %.15lf\n", i, eps);
		i++;
	}

	fclose(EPS);
	fclose(DEGREE2);
	fclose(FACTERR);

	return 0;
}