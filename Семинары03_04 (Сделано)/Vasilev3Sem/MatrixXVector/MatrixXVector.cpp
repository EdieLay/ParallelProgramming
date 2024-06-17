// MatrixXVector.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>

#define SIZE 8000

using namespace std;

double mtrx[SIZE][SIZE];
double result[SIZE]{};
double vectr[SIZE];

int main()
{
    double ctr = 0;
    for (int i = 0; i < SIZE; i++)
    {
        vectr[i] = i;
        for (int j = 0; j < SIZE; j++)
            mtrx[i][j] = ctr++;
    }

    

    for (int k = 6; k < 7; k *= 2)
    {
        omp_set_num_threads(k);
        double time = omp_get_wtime();

#pragma omp parallel for
        for (int i = 0; i < SIZE; i++)
            for (int j = 0; j < SIZE; j++)
                result[i] += mtrx[i][j] * vectr[j];

        time = omp_get_wtime() - time;

        cout << time << endl;

    }
}


