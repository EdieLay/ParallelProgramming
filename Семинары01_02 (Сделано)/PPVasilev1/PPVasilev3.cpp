#define _USE_MATH_DEFINES
#include <iostream>
#include <omp.h>
#include <math.h>


using namespace std;
double pi = 0;
int precision = 11;
double limit = 0.0000000001;
int main()
{
	
	omp_set_num_threads(16);

	double wtime;
	wtime = omp_get_wtime();	
#pragma omp parallel for reduction (+:pi)
	for (int j = 0; j < omp_get_num_threads(); j++)
	{
		for (int i = j * 2 + 1; ; i = i + 2 * omp_get_num_threads())
		{
			double first = 4.0 / (2.0 * i - 1.0);
			double second = 4.0 / (2.0 * i + 1.0);
			if (first < limit)
				break;
			pi += first;
			if (second < limit)
				break;
			pi -= second;
		}
	}
	// ---------
	wtime = omp_get_wtime() - wtime;

	printf("Real PI: %.*g\n", precision, M_PI);
	printf("Comp PI: %.*g\n", precision, pi);
	cout << "Elapsed time: " << wtime << endl;
}