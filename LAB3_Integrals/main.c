#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>

#define GOLDENRATIO 1.61803398875

//Функция
double f(double x)
{
	return x * x * x * x - 3.2 * x * x * x + 2.5 * x * x - 7 * x - 1.5;
}

//Первообразная
double F(double x)
{
    return 0.2 * x * x * x * x * x - 0.8 * x * x * x * x + ((5.0) / (6.0)) * x * x * x - 3.5 * x * x - 1.5 * x;
}

double NewtonLeibniz(double a, double b)
{
	return F(b) - F(a);
}

//Вторая производная от функции f
double d2f(double x)
{
	return 12 * x * x - 19.2 * x + 5;
}

double d2f_max(double a, double b)
{
	double accuracy=0.00001; 
	double x1, x2; 
	while (fabs(b - a) > accuracy) {
		x1 = b - (b - a) / GOLDENRATIO;
		x2 = a + (b - a) / GOLDENRATIO;
		if (d2f(x1) <= d2f(x2))
			a = x1;
		else
			b = x2;
	} 

	printf("MAXMIMUM OF f''(x): %lf\n", d2f((a + b) / 2));

	return d2f((a + b) / 2);
}

//Вычисление теоретической погрешности
double methodTheoreticalError(double h, double a, double b)
{
	return d2f_max(a, b) * ((b - a) / 24) * h * h;
}


void MiddleRectangularMethod(double a, double b, double eps)
{
	double errorFact = 0, errorTheor = 0;
	double h;

	int n = 10;

	double err;
	int runge = 2;
	double x = a;
	double integralCurr = 0;
	double integralPrev = 0;
	double integral = 0;

	do {
		integralPrev = integral;
		integral = 0;
		h = (b - a) / ((double)n);
		for (int i = 0; i < n; i++)
		{
			x += h;
			integral += f(x - h / 2.0) * h;
		}
		printf("F = %lf\n", NewtonLeibniz(a, b));
		printf("integral = %lf\n", integral);
		err = integral - (NewtonLeibniz(a, b));
		printf("err = %.15lf\n", err);
		printf("iter count = %d\n", n);
		n *= runge;
		integralCurr = integral;
		x = a;
	} while (fabs(integralPrev - integralCurr) / 3.0 > eps);

	errorFact = fabs(NewtonLeibniz(a, b) - integralCurr);
	errorTheor = methodTheoreticalError(h, a, b);

	printf("FACT error: %.15lf\n", errorFact);
	printf("THEOR error: %.15lf", errorTheor);
}

int main()
{
	double a, b, h, eps;

	printf("Enter [a,b]: \n");
	scanf("%lf%lf", &a, &b);
	printf("Enter epsilon for Runge Rule:\n");
	scanf("%lf", &eps);
	
	printf("--------------------------------------------\n");

	MiddleRectangularMethod(a, b, eps);

	return 0;
}

