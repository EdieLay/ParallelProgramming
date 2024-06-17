#include <windows.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

#define N 110000000

double a[N + 1], b[N + 1];

int i;
double gDotProduct = 0;

int main()
{
	// initialize vectors 
	for (i = 0; i < N; i++)
	{
		a[i] = 1.034;
		b[i] = 1.057;
	}

	omp_set_num_threads(5);

	double wtime;
	wtime = omp_get_wtime();
	// код, который нужно распараллелить
#pragma omp parallel for
	for (i = 0; i < N; i++)
	{
		gDotProduct += a[i] * b[i];
	}
	// ---------
	wtime = omp_get_wtime() - wtime;

	printf("Computed value of vector sum: ");
	//print dot product  
	printf("sum = %f\n", gDotProduct);
	//print delta time
	printf("time = %g\n", wtime);
}
