#include "Massive.h"
#include <algorithm>
#include <fstream>
#include "mpi.h"


using namespace std;


void Massive::fillingMAssive(int* array) {
	ifstream in;
	cout << "Process filling massive is started.\n";
	string fileName = config.getDataFileName();
	in.open(fileName);
	if (in.is_open()) {
		cout << "File \"" << fileName << "\" is open\n";
		cout << "Filling array...\n";
#pragma omp parallel shared(array) num_threads(config.getNumberOfThreads())
		{
#pragma omp for
			for (long i = 0; i < config.getSizeOfArray(); i++) {
				in >> array[i];
			}
		}
		cout << "Array is filling.\n";
	}
	else {
		cout << "File \"" << fileName << "\" is not open. Check if this file exists.\n";
	}
	cout << "Process filling massive is ended.\n";
}

unsigned long Massive::sumOfArray(int* array) {
	cout << "Process solve a sum of Array is started...\n";
	unsigned long sum =0;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		sum += array[i];
	}
	end_sync = clock();
	seconds = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Process solve a sum of Array is ended.\n";
	cout << "Execution time is " << seconds << " seconds.\n" << endl;
	return sum;
}

unsigned long Massive::sumOfArrayAsyncMPI(int* array) {
	cout << "Process solve a sum of Array is started...\n";
	unsigned long sum = 0;
	unsigned long localSum = 0;

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Разделение массива между процессами
	unsigned long localSize = config.getSizeOfArray() / size;
	unsigned long startIndex = localSize * rank;
	unsigned long endIndex = startIndex + localSize;

	// Суммирование локальных частей массива
	for (unsigned long i = startIndex; i < endIndex; i++) {
		localSum += array[i];
	}

	// Обмен и суммирование результатов
	MPI_Allreduce(&localSum, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

	if (rank == 0) {
		// Вывод результатов только из процесса с рангом 0
		end_sync = clock();
		seconds = (end_sync - start_sync) / CLOCKS_PER_SEC;
		cout << "Process solve a sum of Array is ended.\n";
		cout << "Execution time is " << seconds << " seconds.\n" << endl;
	}

	return sum;
}


unsigned long Massive::sumOfArrayAsync(int* array) {
	cout << "Process solve a sum of Array is started...\n";
	unsigned long sum = 0;
	start_sync = clock();
#pragma omp parallel shared(array) num_threads(threads)
	{
#pragma omp for
		for (long i = 0; i < config.getSizeOfArray(); i++) {
			sum += array[i];
		}
	}
	end_sync = clock();
	seconds = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Process solve a sum of Array is ended.\n";
	cout << "Execution time is " << seconds << " seconds.\n" << endl;
	return sum;
}

long Massive::partition(int* array, long start, long end)
{
	long pivot = array[end];

	long pIndex = start;

	for (long i = start; i < end; i++)
	{
		if (array[i] <= pivot)
		{
			swap(array[i], array[pIndex]);
			pIndex++;
		}
	}

	swap(array[pIndex], array[end]);

	return pIndex;
}

long Massive::partitionAsync(int* array, long start, long end)
{
	long pivot = array[end];

	long pIndex = start;
#pragma omp parallel shared (array) num_threads(threads)
	{
#pragma omp for
		for (long i = start; i < end; i++)
		{
			if (array[i] <= pivot)
			{
				swap(array[i], array[pIndex]);
				pIndex++;
			}
		}

	}
	swap(array[pIndex], array[end]);
	return pIndex;
}

void Massive::quicksort(int* array, long start, long end)
{
	long stackSize = end - start + 1;
	long* stack = new long[stackSize];
	long top = -1;

	stack[++top] = start;
	stack[++top] = end;

	while (top >= 0) {
		end = stack[top--];
		start = stack[top--];

		if (start >= end) {
			continue;
		}

		long pivot = partition(array, start, end);

		if (pivot - 1 > start) {
			stack[++top] = start;
			stack[++top] = pivot - 1;
		}

		if (pivot + 1 < end) {
			stack[++top] = pivot + 1;
			stack[++top] = end;
		}
	}

	delete[] stack;
}

void Massive::quicksortAsync(int* array, long start, long end)
{
	long stackSize = end - start + 1;
	long* stack = new long[stackSize];
	long top = -1;
	stack[++top] = start;
	stack[++top] = end;
#pragma parallel shared (stack,array,start,end) num_threads(threads)
	{
#pragma omp while
		while (top >= 0) {
			end = stack[top--];
			start = stack[top--];

			if (start >= end) {
				continue;
			}

			long pivot = partitionAsync(array, start, end);

			if (pivot - 1 > start) {
				stack[++top] = start;
				stack[++top] = pivot - 1;
			}

			if (pivot + 1 < end) {
				stack[++top] = pivot + 1;
				stack[++top] = end;
			}
		}
	}
	delete[] stack;
}