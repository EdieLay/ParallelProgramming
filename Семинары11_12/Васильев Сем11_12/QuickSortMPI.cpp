#include <iostream>
#include <stdio.h> 
#include "mpi.h" 
#include <random>

using namespace std;

void ParallelHyperQuickSort(double* pProcData, int ProcDataSize);
void PivotDistribution(double* pProcData, int ProcDataSize, int Dim, int Mask, int Iter, double* pPivot);
void ProcessInitialization(double*& pProcData, int& ProcDataSize);
void ProcessTermination(double* pProcData, int ProcDataSize);
void LocalDataSort(double* pProcData, int ProcDataSize);
void QuickSort(double* A, int i1, int i2);
int GetProcDataDivisionPos(double* pProcData, int ProcDataSize, int Pivot);
void DataMerge(double* pMergeData, double* pData, int DataSize, double* pRecvData, int RecvDataSize);

int ProcRank; // Rank of current process 
int ProcNum; // Number of processes 
int main(int argc, char* argv[]) {
	double* pProcData; // Блок данных для процесса 
	int ProcDataSize; // Размер блока данных 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	// Инициализация данных и их распределение между процессами 
	ProcessInitialization(pProcData, ProcDataSize);
	double wtime = 0;
	if (ProcRank == 0)
	{
		wtime = MPI_Wtime();
	}
	// параллельная сортировка 
	ParallelHyperQuickSort(pProcData, ProcDataSize);
	if (ProcRank == 0)
	{
		wtime = MPI_Wtime() - wtime;
		cout << "wtime = " << wtime << endl;
	}
	// Завершение вычислений процесса 
	ProcessTermination(pProcData, ProcDataSize);
	MPI_Finalize();
}

void ProcessInitialization(double*& pProcData, int& ProcDataSize)
{
	ProcDataSize = 1000000 / ProcNum;
	pProcData = new double[ProcDataSize];
	random_device crypto_random_generator;
	uniform_int_distribution<int> int_distribution(100, 999);

	for (int i = 0; i < ProcDataSize; i++) {
		int result = int_distribution(crypto_random_generator);
		pProcData[i] = result;
	}
}

void ProcessTermination(double* pProcData, int ProcDataSize)
{
	// крашится при delete[] pProcData;
}

void ParallelHyperQuickSort(double* pProcData, int ProcDataSize) {
	MPI_Status status;
	int CommProcRank; // ранг процессора, с которым выполняется взаимодействие 
	double* pData, // часть блока, остающаяся на процессоре 
		* pSendData, // часть блока, передаваемая процессору CommProcRank 
		* pRecvData, // часть блока, получаемая от процессора CommProcRank 
		* pMergeData; // блок данных, получаемый после слияния 
	int DataSize, SendDataSize, RecvDataSize, MergeDataSize;
	int HypercubeDim = (int)ceil(log(ProcNum) / log(2)); // размерность гиперкуба 
	int Mask = ProcNum;
	double Pivot;
	// первоначальная сортировка блоков данных на каждом процессоре 
	LocalDataSort(pProcData, ProcDataSize);
	// итерации обобщенной быстрой сортировки 
	for (int i = HypercubeDim; i > 0; i--) {
		// определение ведущего значения и его рассылка всем процессорам 
		PivotDistribution(pProcData, ProcDataSize, HypercubeDim, Mask, i, &Pivot);
		Mask = Mask >> 1;
		// определение границы разделения блока 
		int pos = GetProcDataDivisionPos(pProcData, ProcDataSize, Pivot);
		// разделение блока на части 
		if (((ProcRank & Mask) >> (i - 1)) == 0) { // старший бит = 0 
			pSendData = &pProcData[pos + 1];
			SendDataSize = ProcDataSize - pos - 1;
			if (SendDataSize < 0) SendDataSize = 0;
			CommProcRank = ProcRank + Mask;
			pData = &pProcData[0];
			DataSize = pos + 1;
		}
		else { // старший бит = 1 
			pSendData = &pProcData[0];
			SendDataSize = pos + 1;
			if (SendDataSize > ProcDataSize) SendDataSize = pos;
			CommProcRank = ProcRank - Mask;
			pData = &pProcData[pos + 1];
			DataSize = ProcDataSize - pos - 1;
			if (DataSize < 0) DataSize = 0;
		}
		// пересылка размеров частей блоков данных 
		MPI_Sendrecv(&SendDataSize, 1, MPI_INT, CommProcRank, 0,
			&RecvDataSize, 1, MPI_INT, CommProcRank, 0, MPI_COMM_WORLD, &status);
		// пересылка частей блоков данных 
		pRecvData = new double[RecvDataSize];
		MPI_Sendrecv(pSendData, SendDataSize, MPI_DOUBLE,
			CommProcRank, 0, pRecvData, RecvDataSize, MPI_DOUBLE,
			CommProcRank, 0, MPI_COMM_WORLD, &status);

		// слияние частей 
		MergeDataSize = DataSize + RecvDataSize;
		pMergeData = new double[MergeDataSize];
		DataMerge(pMergeData, pData, DataSize, pRecvData, RecvDataSize);
		delete[] pProcData;
		delete[] pRecvData;
		pProcData = pMergeData;
		ProcDataSize = MergeDataSize;
	}
}



void LocalDataSort(double* pProcData, int ProcDataSize)
{
	QuickSort(pProcData, 0, ProcDataSize);
}

void QuickSort(double* A, int i1, int i2) {
	if (i1 < i2) {
		double pivot = A[i1];
		int is = i1;
		for (int i = i1 + 1; i < i2; i++)
			if (A[i] <= pivot) {
				is = is + 1;
				swap(A[is], A[i]);
			}
		swap(A[i1], A[is]);
		QuickSort(A, i1, is);
		QuickSort(A, is + 1, i2);
	}
}

int GetProcDataDivisionPos(double* pProcData, int ProcDataSize, int Pivot)
{
	for (int i = 0; i < ProcDataSize; i++)
	{
		if (pProcData[i] > Pivot) return i;
	}
}

void DataMerge(double* pMergeData, double* pData, int DataSize, double* pRecvData, int RecvDataSize)
{
	int i = 0, j = 0, k = 0; 
    while (i < DataSize && j < RecvDataSize) 
    { 
        if (pData[i] < pRecvData[j]) 
            pMergeData[k++] = pData[i++]; 
        else
            pMergeData[k++] = pRecvData[j++]; 
    } 
  
    while (i < DataSize) 
        pMergeData[k++] = pData[i++]; 
  
    while (j < RecvDataSize) 
        pMergeData[k++] = pRecvData[j++]; 
}

void PivotDistribution(double* pProcData, int ProcDataSize, int Dim, int Mask, int Iter, double* pPivot)
{
	MPI_Group WorldGroup;
	MPI_Group SubcubeGroup; // группа процессов - подгиперкуб 
	MPI_Comm SubcubeComm; // коммуникатор подгиперкуба 
	int j = 0;
	int GroupNum = ProcNum / (int)pow(2, Dim - Iter);
	int* ProcRanks = new int[GroupNum];
	// формирование списка рангов процессов для гиперкуба 
	int StartProc = ProcRank - GroupNum;
	if (StartProc < 0) 
		StartProc = 0;

	int EndProc = ProcRank + GroupNum;
	if (EndProc > ProcNum) 
		EndProc = ProcNum;

	for (int proc = StartProc; proc < EndProc; proc++) 
	{
		if ((ProcRank & Mask) >> (Iter) == (proc & Mask) >> (Iter)) 
		{
			ProcRanks[j++] = proc;
		}
	}
	//объединение процессов подгиперкуба в одну группу 
	MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
	MPI_Group_incl(WorldGroup, GroupNum, ProcRanks, &SubcubeGroup);
	MPI_Comm_create(MPI_COMM_WORLD, SubcubeGroup, &SubcubeComm);
	// поиск и рассылка ведущего элемента всем процессам подгиперкуба 
	if (ProcRank == ProcRanks[0])
		*pPivot = pProcData[(ProcDataSize) / 2];

	MPI_Bcast(pPivot, 1, MPI_DOUBLE, 0, SubcubeComm);
	MPI_Group_free(&SubcubeGroup);
	MPI_Comm_free(&SubcubeComm);
	delete[] ProcRanks;
}