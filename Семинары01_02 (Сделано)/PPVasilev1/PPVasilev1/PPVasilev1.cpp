#include <iostream>
#include <clocale>
#include <omp.h>
#include <cmath>

using namespace std;

bool is_prime(int num)
{
    int num_sqrt = sqrt(num);
    for (int i = 2; i <= num_sqrt; i++)
    {
        if (num % i == 0)
            return false;
    }
    return true;
}

int main()
{
    int tid;
    setlocale(LC_ALL, "Russian");

#ifdef _OPENMP 
    cout << "Поддержка OpenMP включена" << endl;
#else
    cout << "OpenMP не поддерживается" << endl;
#endif

    omp_set_num_threads(8);

    int sum = 0;

    double wtime;
    wtime = omp_get_wtime();

#pragma omp parallel for reduction (+:sum)
        for (int i = 2; i < 100000; i++)
            if (is_prime(i))
                sum += i;

    wtime = omp_get_wtime() - wtime;

    cout << sum << endl;
    cout << wtime << endl;
}
