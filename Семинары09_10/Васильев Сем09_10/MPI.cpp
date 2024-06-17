#include <iostream>
#include "mpi.h"

using namespace std;

int ProcNum; // количество процессоров
int ProcRank; // номер процессора, на котором выполняются описываемые действия
int PrevProc, NextProc; // номера соседних процессоров, содержащих предшествующую и следующую полосы
int M; // количество строк в полосе (без учета продублированных граничных строк)
int N; // количество внутренних узлов в строке сетки (т.е. всего в строке N+2 узла)
double** f;
double** u;

int main(int argc, char* argv[])
{
	double dmax, eps = 0.0001;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	PrevProc = ProcRank - 1;
	NextProc = ProcRank + 1;
	N = 3000;
	M = N / ProcNum;
	double h = 1.0 / ((double)N + 1.0);
	f = new double* [N + 2];
	u = new double* [M + 2];
	for (int i = 0; i < N + 2; i++)
	{
		f[i] = new double[N + 2];
		for (int j = 0; j < N + 2; j++)
		{
			f[i][j] = 1;
		}
	}
	for (int i = 0; i < M + 2; i++)
	{
		u[i] = new double[N + 2];
		for (int j = 0; j < N + 2; j++)
		{
			u[i][j] = (i + M * ProcRank) * N + j;
		}
	}

	double wtime = MPI_Wtime();
	do {
		MPI_Status status;
		// обмен граничных строк полос с соседями 
		if (ProcRank > 0 && ProcRank < ProcNum - 1)
		{
			MPI_Sendrecv(u[M], N + 2, MPI_DOUBLE, NextProc, 1,
				u[0], N + 2, MPI_DOUBLE, PrevProc, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
			MPI_Sendrecv(u[1], N + 2, MPI_DOUBLE, PrevProc, 1,
				u[M + 1], N + 2, MPI_DOUBLE, NextProc, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
		}
		else if (ProcRank < ProcNum - 1)
		{
			MPI_Send(u[M], N + 2, MPI_DOUBLE, NextProc, 1, MPI_COMM_WORLD);
			MPI_Recv(u[M + 1], N + 2, MPI_DOUBLE, NextProc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
		else if (ProcRank > 0)
		{
			MPI_Send(u[1], N + 2, MPI_DOUBLE, PrevProc, 1, MPI_COMM_WORLD);
			MPI_Recv(u[0], N + 2, MPI_DOUBLE, PrevProc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}

		double temp, d, dm;
		int i;
		dmax = 0; // максимальное изменение значений u  
		dm = 0;
		for (i = 1; i < M + 1; i++) {
			for (int j = 1; j < N + 1; j++) {
				temp = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] +
					u[i][j - 1] + u[i][j + 1] - h * h * f[i][j]);
				d = fabs(temp - u[i][j]);
				if (dm < d) dm = d;
			}
		}

		// Поиск максимального dmax среди dm на процессе ранга 0
		MPI_Reduce(&dm, &dmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		// рассылка dmax по всем процессам с процесса ранга 0
		MPI_Bcast(&dmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} while (dmax > eps); // eps - точность решения 
	wtime = MPI_Wtime() - wtime;
	if (ProcRank == 0)
		cout << wtime << endl;

	MPI_Finalize();
}
