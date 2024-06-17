
#define _USE_MATH_DEFINES
#include <iostream>
#include <omp.h>
#include <math.h>


using namespace std;

int main()
{
	double pi = 0;
	omp_set_num_threads(5);

	double wtime;
	wtime = omp_get_wtime();

	int precision = 11;
	
#pragma omp parallel for reduction (+:pi)
	for (int i = 1; 1.0 / (2.0*(double)i - 1.0) > 1.0 / pow(10.0, (double)precision + 3.0); i = i+2)
	{
		pi += 1.0 / (2.0 * i - 1.0) - 1.0 / (2.0 * i + 1.0);
	}
	pi *= 4;
	// ---------
	wtime = omp_get_wtime() - wtime;

	printf("Real PI: %.*g\n", precision, M_PI);
	printf("Comp PI: %.*g\n", precision, pi);
	cout << "Elapsed time: " << wtime << endl;
}
