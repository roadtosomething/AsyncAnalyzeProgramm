#include "Operation.h"
#include "omp.h"
#include "Config.h"
#include <iostream>
#include <time.h>

using namespace std;
Config config;
int threads = config.getNumberOfThreads();
clock_t start_sync;
clock_t end_sync;

void Operation::addition(int* array_a, int* array_b) {
	cout << "Synchronous operation addition is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] + array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation addition is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::addition(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation addition is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] + matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation addition is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void subtraction(int* array_a, int* array_b) {
	cout << "Synchronous operation subtraction is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] - array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation subtraction is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void subtraction(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation subtraction is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] - matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation subtraction is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void multiplication(int* array_a, int* array_b) {
	cout << "Synchronous operation multiplication is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] * array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation multiplication is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void multiplication(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation multiplication is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] * matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation multiplication is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds";
}

void division(int* array_a, int* array_b) {
	cout << "Synchronous operation division is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] / array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation division is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void division(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation division is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] / matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation division is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::additionAsync(int* array_a, int* array_b) {
	cout << "Asynchronous operation addition is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (long i = 0; i < config.getSizeOfArray(); i++) {
			array_a[i] + array_b[i];
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation addition is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void additionAsync(int** matrix_a, int** matrix_b) {
	cout << "Asynchronous operation addition is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (int j = 0; j < config.getSizeOfMatrix(); j++) {
#pragma omp for
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				matrix_a[j][i] + matrix_b[j][i];
			}
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation addition is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void subtractionAsync(int* array_a, int* array_b) {
	cout << "Asynchronous operation subtraction is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (long i = 0; i < config.getSizeOfArray(); i++) {
			array_a[i] - array_b[i];
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation subtraction is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void subtractionAsync(int** matrix_a, int** matrix_b) {
	cout << "Asynchronous operation subtraction is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (int j = 0; j < config.getSizeOfMatrix(); j++) {
#pragma omp for
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				matrix_a[j][i] - matrix_b[j][i];
			}
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation subtraction is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void multiplicationAsync(int* array_a, int* array_b) {
	cout << "Asynchronous operation multiplication is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (long i = 0; i < config.getSizeOfArray(); i++) {
			array_a[i] * array_b[i];
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation multiplication is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void multiplicationAsync(int** matrix_a, int** matrix_b) {
	cout << "Asynchronous operation multiplication is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (int j = 0; j < config.getSizeOfMatrix(); j++) {
#pragma omp for
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				matrix_a[j][i] * matrix_b[j][i];
			}
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation multiplication is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void divisionAsync(int* array_a, int* array_b) {
	cout << "Asynchronous operation division is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (long i = 0; i < config.getSizeOfArray(); i++) {
			array_a[i] / array_b[i];
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation division is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void divisionAsync(int** matrix_a, int** matrix_b) {
	cout << "Asynchronous operation division is started" << endl;
	start_sync = clock();
#pragma omp parallel shared(array_a, array_b) num_threads(threads)
	{
#pragma omp for
		for (int j = 0; j < config.getSizeOfMatrix(); j++) {
#pragma omp for
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				matrix_a[j][i] / matrix_b[j][i];
			}
		}
	}
	end_sync = clock();
	cout << "Asynchronous operation division is ended" << endl;
	double time = (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}