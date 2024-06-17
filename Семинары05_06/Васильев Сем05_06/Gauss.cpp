#include <iostream>
#include <random>
#include <omp.h>

#define SIZE 3000

using namespace std;

double* SolveGauss(double** a, double* b, int n);
void FindMaxAndSwap(double** a, double* b, int k, int n);
void FillRandomMatrix(double** a, double* b, int n);

int main()
{
    double** a = new double* [SIZE];
    for (int i = 0; i < SIZE; i++)
        a[i] = new double[SIZE];
    double* b = new double[SIZE];
    FillRandomMatrix(a, b, SIZE);

    double time = omp_get_wtime();
    double* x = SolveGauss(a, b, SIZE);
    time = omp_get_wtime() - time;
    cout << time << endl;
}


double* SolveGauss(double** a, double* b, int n)
{
    omp_set_num_threads(4);
    for (int k = 0; k < n; k++) // прямой ход (по столбцам)
    {
        FindMaxAndSwap(a, b, k, n);

#pragma omp parallel for
        for (int j = k + 1; j < n; j++) // идем по строкам
        {
            double d = a[j][k] / a[k][k]; // для текущей строки вычисляем делитель для приведения к верхнетреугольному виду
            for (int i = k; i < n; i++)
            {
                a[j][i] = a[j][i] - d * a[k][i]; // вычитаем текущую строку, умноженную на делитель, из всех строк ниже её
            }
            b[j] = b[j] - d * b[k]; // для столбца значений то же самое
        }
    }

    double* x = new double[n + 1]; // массив, где будет храниться решение

#pragma omp parallel for
    for (int k = n-1; k >= 0; k--) // обратный ход
    {
        double d = 0;
        for (int j = k + 1; j <= n; j++)
        {
            double s = a[k][j] * x[j];
            d = d + s;
        }
        x[k] = (b[k] - d) / a[k][k];
    }

    return x;
}

void FindMaxAndSwap(double** a, double* b, int k, int n)
{
    double max = abs(a[k][k]);
    int max_ind = k;
    for (int i = k + 1; i < n; i++)
    {
        double cur = abs(a[i][k]);
        if (cur > max)
        {
            max = cur;
            max_ind = i;
        }
    }
    if (max_ind != k)
    {
        double temp;
        for (int i = k; i < n; i++)
        {
            temp = a[k][i];
            a[k][i] = a[max_ind][i];
            a[max_ind][i] = temp;
        }
        temp = b[k];
        b[k] = b[max_ind];
        b[max_ind] = temp;
    }
}

void FillRandomMatrix(double** a, double* b, int n)
{
    random_device crypto_random_generator;
    uniform_int_distribution<int> int_distribution(-1000, 1000);
    int result;

    for (int i = 0; i < n; i++) {
        for (int j = 0; i < n; i++)
        {
            result = int_distribution(crypto_random_generator);
            a[i][j] = result;
        }
        result = int_distribution(crypto_random_generator);
        b[i] = result;
    }
}