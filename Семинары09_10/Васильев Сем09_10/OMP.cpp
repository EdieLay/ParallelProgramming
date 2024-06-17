#include <iostream>
#include <omp.h>
#define N_NUM 5
#define THREAD_NUM_NUM 4

using namespace std;

int main()
{
	int NValues[N_NUM] = { 10, 100, 1000, 2000, 3000 };
	int threadNum[THREAD_NUM_NUM] = { 1, 2, 4, 6 };
	int N;
	double times[N_NUM][THREAD_NUM_NUM];
	for (int k = 0; k < N_NUM; k++)
	{
		N = NValues[k];
		double** f = new double* [N + 2];
		double** u = new double* [N + 2];
		for (int i = 0; i < N + 2; i++)
		{
			f[i] = new double[N + 2];
			u[i] = new double[N + 2];
			for (int j = 0; j < N + 2; j++)
			{
				f[i][j] = 1;
				u[i][j] = i * N + j;
			}
		}
		omp_lock_t dmax_lock;
		omp_init_lock(&dmax_lock);
		double dmax, eps = 0.0001, h = 1.0 / ((double)N + 1.0);

		for (int l = 0; l < THREAD_NUM_NUM; l++)
		{
			omp_set_num_threads(threadNum[l]);
			double wtime = omp_get_wtime();
			do {
				double temp, d, dm;
				int i;
				dmax = 0; // максимальное изменение значений u 
#pragma omp parallel for shared(u,N,dmax)private(i,temp,d,dm) 
				for (i = 1; i < N + 1; i++) {
					dm = 0;
					for (int j = 1; j < N + 1; j++) {
						temp = u[i][j];
						u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] +
							u[i][j - 1] + u[i][j + 1] - h * h * f[i][j]);
						d = fabs(temp - u[i][j]);
						if (dm < d) dm = d;
					}
					omp_set_lock(&dmax_lock);
					if (dmax < dm) dmax = dm;
					omp_unset_lock(&dmax_lock);
				} // конец параллельной области 
			} while (dmax > eps);
			wtime = omp_get_wtime() - wtime;
			times[k][l] = wtime;
		}
	}

	for (int j = 0; j < THREAD_NUM_NUM; j++)
	{
		cout << endl << threadNum[j] << " Threads" << endl;
		for (int i = 0; i < N_NUM; i++)
		{
			cout << times[i][j] << endl;
		}
	}
}

