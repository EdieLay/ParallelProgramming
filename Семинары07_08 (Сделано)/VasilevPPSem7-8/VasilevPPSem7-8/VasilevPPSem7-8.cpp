// VasilevPPSem7-8.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <stdio.h> 
#include "mpi.h" 
#define SIZE 120000

double vector[SIZE];

int main(int argc, char* argv[]) {
	for (int i = 0; i < SIZE; i++)
	{
		vector[i] = i;
	}
	int ProcNum, ProcRank, RecvRank;
	double sum = 0, procSum = 0;
	double wtime = 0;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int k = SIZE / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);

	if (ProcRank == 0)
		wtime = MPI_Wtime();

	MPI_Bcast(vector, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = i1; i < i2; i++)
		procSum += vector[i];
	MPI_Reduce(&procSum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		wtime = MPI_Wtime() - wtime;
		std::cout << sum << std::endl;
		std::cout << wtime << std::endl;
	}
	MPI_Finalize();
	return 0;
}


