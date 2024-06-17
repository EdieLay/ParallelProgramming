#include <iostream>
#include <stdio.h> 
#include "mpi.h" 

#define SIZE 1000

using namespace std;

double mtrx[SIZE][SIZE];
double result[SIZE]{};
double procResult[SIZE]{};
double vectr[SIZE];

int main(int argc, char* argv[])
{
    double ctr = 0;
    for (int i = 0; i < SIZE; i++)
    {
        vectr[i] = i;
        for (int j = 0; j < SIZE; j++)
            mtrx[i][j] = ctr++;
    }
    int ProcNum, ProcRank, RecvRank;
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

    MPI_Bcast(vectr, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(mtrx, SIZE*SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = i1; i < i2; i++)
        for (int j = 0; j < SIZE; j++)
            procResult[i] += mtrx[i][j] * vectr[j];

    MPI_Reduce(procResult, result, SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (ProcRank == 0)
    {
        wtime = MPI_Wtime() - wtime;
        cout << wtime << endl;
        //for (int i = 0; i < SIZE; i++)
        //    cout << result[i] << " ";
    }

    MPI_Finalize();
    return 0;
}

