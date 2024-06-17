// MatrixXVector.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>

#define SIZE 4000

using namespace std;

double mtrx_a[SIZE][SIZE];
double mtrx_b[SIZE][SIZE];
double result[SIZE][SIZE]{};

int main()
{
    double ctr = 0;
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            mtrx_a[i][j] = ctr++;
            mtrx_b[i][j] = SIZE - ctr;
        }
    }



    for (int num = 1; num < 7; num += 2)
    {
        omp_set_num_threads(num);
        
        for (int BLOCK_SIZE = 50; BLOCK_SIZE < 170; BLOCK_SIZE += 50)
        {
            int GRID_SIZE = SIZE / BLOCK_SIZE;
            double time = omp_get_wtime();

#pragma omp parallel for
            for (int n = 0; n < GRID_SIZE; n++)
                for (int m = 0; m < GRID_SIZE; m++)
                    for (int iter = 0; iter < GRID_SIZE; iter++)
                        for (int i = n * BLOCK_SIZE; i < (n + 1) * BLOCK_SIZE; i++)
                            for (int j = m * BLOCK_SIZE; j < (m + 1) * BLOCK_SIZE; j++)
                                for (int k = iter * BLOCK_SIZE; k < (iter + 1) * BLOCK_SIZE; k++)
                                    result[i][j] += mtrx_a[i][k] * mtrx_b[k][j];

            time = omp_get_wtime() - time;

            cout << time << "\t";
        }
        cout << endl;
        if (1 == num) num--;
    }
}


