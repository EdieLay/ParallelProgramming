#include <iostream>
#include <random>
#include <omp.h>

using namespace std;

pair<int, int> partition(int* arr, int l, int r, int pivot);
void quick_sort_rec(int* arr, int l, int r);
void quick_sort(int* arr, int len);
void FillRandomArray(int* arr, int len);


int main()
{
    int threads_num[4] = { 1, 2, 4, 8 };
    int len[9] = { 1000, 2000, 5000, 10000, 20000, 50000, 100000, 500000, 1000000 };

    for (int i = 0; i < 4; i++)
    {
        omp_set_num_threads(threads_num[i]);
        cout << threads_num[i] << " threads" << endl;
        for (int j = 0; j < 9; j++)
        {
            int* arr = new int[len[j]];
            FillRandomArray(arr, len[j]);
            double wtime = omp_get_wtime();
            quick_sort(arr, len[j]);
            wtime = omp_get_wtime() - wtime;
            cout << wtime << endl;
        }
    }

    return 0;
}

void FillRandomArray(int* arr, int len)
{
    random_device crypto_random_generator;
    uniform_int_distribution<int> int_distribution(100, 999);

    for (int i = 0; i < len; i++) {
        int result = int_distribution(crypto_random_generator);
        arr[i] = result;
    }
}

pair<int, int> partition(int* arr, int l, int r, int pivot) // l - left (включительно), r - right (не включительно)
{
    int e = l; // equal
    int g = l; // greater
    int n = l; // now
    while (n < r)
    {
        if (arr[n] < pivot)
        {
            swap(arr[n], arr[g]); // свапаем текущий с большим
            swap(arr[g], arr[e]); // на месте большего - текущий, поэтому свапаем текущий с равным
            e++;
            g++;
        }
        else if (arr[n] == pivot)
        {
            swap(arr[n], arr[g]);
            g++;
        }
        n++;
    }
    return make_pair(e, g);
}

void quick_sort_rec(int* arr, int l, int r)
{
    if (l >= r)
        return;
    int pivot = arr[rand() % (r - l) + l];
    pair<int, int> res = partition(arr, l, r, pivot);
    int less = res.first;
    int bigger = res.second;
#pragma omp parallel sections
    {
#pragma omp section
        {
            quick_sort_rec(arr, l, less);
        }
#pragma omp section
        {
            quick_sort_rec(arr, bigger, r);
        }
    }
}

void quick_sort(int* arr, int len)
{
    if (len <= 1)
        return;
    int pivot = arr[rand() % len];
    pair<int, int> res = partition(arr, 0, len, pivot);
    int less = res.first;
    int bigger = res.second;
#pragma omp parallel sections
    {
#pragma omp section
        {
            quick_sort_rec(arr, 0, less);
        }
#pragma omp section
        {
            quick_sort_rec(arr, bigger, len);
        }
    }
}
