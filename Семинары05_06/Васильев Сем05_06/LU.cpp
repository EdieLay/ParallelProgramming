#include <iostream>
#include <clocale>
#include <omp.h>

using namespace std;

#define SIZE 3
double A[SIZE][SIZE];
double U[SIZE][SIZE];
double L[SIZE][SIZE];
double C[SIZE][SIZE];
double b[SIZE];
double x[SIZE];
double y[SIZE];

int main()
{
    A[0][0] = 2; A[0][1] = -1; A[0][2] = 1;
    A[1][0] = 3; A[1][1] = -2; A[1][2] = -3;
    A[2][0] = 1; A[2][1] = 1; A[2][2] = 1;

    b[0] = 0; b[1] = 11; b[2] = 16;

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;

            if (i == j)
                U[i][j] = 1.0;

            C[i][j] = 0.0;

            y[i] = b[i];
            x[i] = b[i];
        }
    }

    omp_set_num_threads(4);

#pragma omp parallel for
    // находим первый столбец L и первую строку U
    for (int i = 0; i < SIZE; i++) {
        L[i][0] = A[i][0];
        U[0][i] = A[0][i] / L[0][0];
    }

    double sum = 0;

    for (int i = 1; i < SIZE; i++) {
#pragma omp parallel for
        for (int j = 1; j < SIZE; j++) {
            if (i >= j) { // нижний треугольник
                sum = 0;
                for (int k = 0; k < j; k++)
                    sum += L[i][k] * U[k][j];

                L[i][j] = A[i][j] - sum;
            }
            else { // верхний
                sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[i][k] * U[k][j];

                U[i][j] = (A[i][j] - sum) / L[i][i];
            }
        }
    }

    cout << "L\n";
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
            cout << L[i][j] << "\t";
        cout << "\n";
    }

    cout << "\nU\n";
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
            cout << U[i][j] << "\t";
        cout << "\n";
    }


    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
#pragma omp parallel for
            for (int k = 0; k < SIZE; k++)
                C[i][j] += L[i][k] * U[k][j];


    cout << "\nC\n";
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
            cout << C[i][j] << "\t";
        cout << "\n";
    }

    for (int i = 0; i < SIZE; i++) {
#pragma omp parallel for
        for (int j = 0; j < i; j++)
            y[i] -= L[i][j] * y[j];
        y[i] /= L[i][i];
    }
#pragma omp parallel for
    for (int i = SIZE - 1; i >= 0; i--) {
        x[i] = y[i];
    }

    for (int i = SIZE - 1; i >= 0; i--) {
#pragma omp parallel for
        for (int j = i + 1; j < SIZE; j++)
            x[i] -= U[i][j] * x[j];
    }

    cout << "\nx\n";
    for (int i = 0; i < SIZE; i++)
        cout << x[i] << "\t";
    cout << "\n";
}
