#include "Operation.h"
#include <omp.h>
#include "Config.h"
#include <iostream>
#include <time.h>
#include "mpi.h"

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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::subtraction(int* array_a, int* array_b) {
	cout << "Synchronous operation subtraction is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] - array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation subtraction is ended" << endl;
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::subtraction(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation subtraction is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] - matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation subtraction is ended" << endl;
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::multiplication(int* array_a, int* array_b) {
	cout << "Synchronous operation multiplication is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] * array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation multiplication is ended" << endl;
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::multiplication(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation multiplication is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] * matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation multiplication is ended" << endl;
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds";
}

void Operation::division(int* array_a, int* array_b) {
	cout << "Synchronous operation division is started" << endl;
	start_sync = clock();
	for (long i = 0; i < config.getSizeOfArray(); i++) {
		array_a[i] / array_b[i];
	}
	end_sync = clock();
	cout << "Synchronous operation division is ended" << endl;
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::division(int** matrix_a, int** matrix_b) {
	cout << "Synchronous operation division is started" << endl;
	start_sync = clock();
	for (long j = 0; j < config.getSizeOfMatrix(); j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			matrix_a[j][i] / matrix_b[j][i];
		}
	}
	end_sync = clock();
	cout << "Synchronous operation division is ended" << endl;
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

//OMP Operation

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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::additionAsync(int** matrix_a, int** matrix_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::subtractionAsync(int* array_a, int* array_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::subtractionAsync(int** matrix_a, int** matrix_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::multiplicationAsync(int* array_a, int* array_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::multiplicationAsync(int** matrix_a, int** matrix_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::divisionAsync(int* array_a, int* array_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

void Operation::divisionAsync(int** matrix_a, int** matrix_b) {
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
	double time = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << time << " seconds\n" << endl;
}

//MPI Operation

void Operation::additionMPI(int* array_a, int* array_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation addition is started" << endl;
	int chunk_size = config.getSizeOfArray() / size;
	int remainder = config.getSizeOfArray() % size;

	int* local_a = new int[chunk_size];
	int* local_b = new int[chunk_size];
	int* local_c = new int[chunk_size];

	MPI_Scatter(array_a, chunk_size, MPI_INT, local_a, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(array_b, chunk_size, MPI_INT, local_b, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < chunk_size; i++) {
		local_c[i] = local_a[i] + local_b[i];
	}

	int* result = nullptr;
	if (rank == 0) {
		result = new int[config.getSizeOfArray()];
	}

	MPI_Gather(local_c, chunk_size, MPI_INT, result, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long i = config.getSizeOfArray() - remainder; i < config.getSizeOfArray(); i++) {
			result[i] = array_a[i] + array_b[i];
		}

		// Print the result
		cout << "MPI operation addition is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		delete[] result;
	}

	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::additionMPI(int** matrix_a, int** matrix_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation addition is started" << endl;
	long chunk_size = config.getSizeOfMatrix() / size;
	int remainder = config.getSizeOfMatrix() % size;

	int** local_a = new int* [chunk_size];
	int** local_b = new int* [chunk_size];
	int** local_c = new int* [chunk_size];

	for (long i = 0; i < chunk_size; i++) {
		local_a[i] = new int[config.getSizeOfMatrix()];
		local_b[i] = new int[config.getSizeOfMatrix()];
		local_c[i] = new int[config.getSizeOfMatrix()];
	}

	MPI_Scatter(&matrix_a[0][0], chunk_size * config.getSizeOfMatrix(), MPI_INT, &local_a[0][0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(&matrix_b[0][0], chunk_size * config.getSizeOfMatrix(), MPI_INT, &local_b[0][0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	for (long j = 0; j < chunk_size; j++) {
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			local_c[j][i] = local_a[j][i] + local_b[j][i];
		}
	}

	int** result = nullptr;
	if (rank == 0) {
		result = new int* [config.getSizeOfMatrix()];
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			result[i] = new int[config.getSizeOfMatrix()];
		}
	}

	MPI_Gather(local_c[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, result[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long j = config.getSizeOfMatrix() - remainder; j < config.getSizeOfMatrix(); j++) {
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				result[j][i] = matrix_a[j][i] + matrix_b[j][i];
			}
		}

		// Print the result
		cout << "MPI operation addition is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			delete[] result[i];
		}
		delete[] result;
	}

	for (int i = 0; i < chunk_size; i++) {
		delete[] local_a[i];
		delete[] local_b[i];
		delete[] local_c[i];
	}
	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::subtractionMPI(int* array_a, int* array_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation subtraction is started" << endl;
	int chunk_size = config.getSizeOfArray() / size;
	int remainder = config.getSizeOfArray() % size;

	int* local_a = new int[chunk_size];
	int* local_b = new int[chunk_size];
	int* local_c = new int[chunk_size];

	MPI_Scatter(array_a, chunk_size, MPI_INT, local_a, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(array_b, chunk_size, MPI_INT, local_b, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < chunk_size; i++) {
		local_c[i] = local_a[i] - local_b[i];
	}

	int* result = nullptr;
	if (rank == 0) {
		result = new int[config.getSizeOfArray()];
	}

	MPI_Gather(local_c, chunk_size, MPI_INT, result, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long i = config.getSizeOfArray() - remainder; i < config.getSizeOfArray(); i++) {
			result[i] = array_a[i] - array_b[i];
		}

		// Print the result
		cout << "MPI operation subtraction is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		delete[] result;
	}

	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::subtractionMPI(int** matrix_a, int** matrix_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation subtraction is started" << endl;
	int chunk_size = config.getSizeOfMatrix() / size;
	int remainder = config.getSizeOfMatrix() % size;

	int** local_a = new int* [chunk_size];
	int** local_b = new int* [chunk_size];
	int** local_c = new int* [chunk_size];

	for (int i = 0; i < chunk_size; i++) {
		local_a[i] = new int[config.getSizeOfMatrix()];
		local_b[i] = new int[config.getSizeOfMatrix()];
		local_c[i] = new int[config.getSizeOfMatrix()];
	}

	MPI_Scatter(matrix_a[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, local_a[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix_b[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, local_b[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	for (int j = 0; j < chunk_size; j++) {
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			local_c[j][i] = local_a[j][i] - local_b[j][i];
		}
	}

	int** result = nullptr;
	if (rank == 0) {
		result = new int* [config.getSizeOfMatrix()];
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			result[i] = new int[config.getSizeOfMatrix()];
		}
	}

	MPI_Gather(local_c[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, result[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long j = config.getSizeOfMatrix() - remainder; j < config.getSizeOfMatrix(); j++) {
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				result[j][i] = matrix_a[j][i] - matrix_b[j][i];
			}
		}

		// Print the result
		cout << "Synchronous operation subtraction is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			delete[] result[i];
		}
		delete[] result;
	}

	for (int i = 0; i < chunk_size; i++) {
		delete[] local_a[i];
		delete[] local_b[i];
		delete[] local_c[i];
	}
	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::multiplicationMPI(int* array_a, int* array_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation multiplication is started..." << endl;
	int chunk_size = config.getSizeOfArray() / size;
	int remainder = config.getSizeOfArray() % size;

	int* local_a = new int[chunk_size];
	int* local_b = new int[chunk_size];
	int* local_c = new int[chunk_size];

	MPI_Scatter(array_a, chunk_size, MPI_INT, local_a, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(array_b, chunk_size, MPI_INT, local_b, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < chunk_size; i++) {
		local_c[i] = local_a[i] * local_b[i];
	}

	int* result = nullptr;
	if (rank == 0) {
		result = new int[config.getSizeOfArray()];
	}

	MPI_Gather(local_c, chunk_size, MPI_INT, result, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long i = config.getSizeOfArray() - remainder; i < config.getSizeOfArray(); i++) {
			result[i] = array_a[i] * array_b[i];
		}

		// Print the result
		cout << "Synchronous operation multiplication is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		delete[] result;
	}

	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::multiplicationMPI(int** matrix_a, int** matrix_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation multiplication is started..." << endl;
	int chunk_size = config.getSizeOfMatrix() / size;
	int remainder = config.getSizeOfMatrix() % size;

	int** local_a = new int* [chunk_size];
	int** local_b = new int* [chunk_size];
	int** local_c = new int* [chunk_size];

	for (int i = 0; i < chunk_size; i++) {
		local_a[i] = new int[config.getSizeOfMatrix()];
		local_b[i] = new int[config.getSizeOfMatrix()];
		local_c[i] = new int[config.getSizeOfMatrix()];
	}

	MPI_Scatter(matrix_a[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, local_a[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix_b[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, local_b[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	for (int j = 0; j < chunk_size; j++) {
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			local_c[j][i] = local_a[j][i] * local_b[j][i];
		}
	}

	int** result = nullptr;
	if (rank == 0) {
		result = new int* [config.getSizeOfMatrix()];
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			result[i] = new int[config.getSizeOfMatrix()];
		}
	}

	MPI_Gather(local_c[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, result[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long j = config.getSizeOfMatrix() - remainder; j < config.getSizeOfMatrix(); j++) {
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				result[j][i] = matrix_a[j][i] * matrix_b[j][i];
			}
		}

		// Print the result
		cout << "Synchronous operation multipliacation is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			delete[] result[i];
		}
		delete[] result;
	}

	for (int i = 0; i < chunk_size; i++) {
		delete[] local_a[i];
		delete[] local_b[i];
		delete[] local_c[i];
	}
	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::divisionMPI(int* array_a, int* array_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation division is started..." << endl;
	int chunk_size = config.getSizeOfArray() / size;
	int remainder = config.getSizeOfArray() % size;

	int* local_a = new int[chunk_size];
	int* local_b = new int[chunk_size];
	double* local_c = new double[chunk_size];

	MPI_Scatter(array_a, chunk_size, MPI_INT, local_a, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(array_b, chunk_size, MPI_INT, local_b, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < chunk_size; i++) {
		local_c[i] = local_a[i] / local_b[i];
	}

	double* result = nullptr;
	if (rank == 0) {
		result = new double[config.getSizeOfArray()];
	}

	MPI_Gather(local_c, chunk_size, MPI_DOUBLE, result, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long i = config.getSizeOfArray() - remainder; i < config.getSizeOfArray(); i++) {
			result[i] = array_a[i] / array_b[i];
		}

		// Print the result
		cout << "MPI operation division is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		delete[] result;
	}

	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}

void Operation::divisionMPI(int** matrix_a, int** matrix_b) {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	cout << "MPI operation division is started..." << endl;
	int chunk_size = config.getSizeOfMatrix() / size;
	int remainder = config.getSizeOfMatrix() % size;

	int** local_a = new int* [chunk_size];
	int** local_b = new int* [chunk_size];
	double** local_c = new double* [chunk_size];

	for (int i = 0; i < chunk_size; i++) {
		local_a[i] = new int[config.getSizeOfMatrix()];
		local_b[i] = new int[config.getSizeOfMatrix()];
		local_c[i] = new double[config.getSizeOfMatrix()];
	}

	MPI_Scatter(matrix_a[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, local_a[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix_b[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, local_b[0], chunk_size * config.getSizeOfMatrix(), MPI_INT, 0, MPI_COMM_WORLD);

	for (int j = 0; j < chunk_size; j++) {
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			local_c[j][i] = local_a[j][i] / local_b[j][i];
		}
	}

	double** result = nullptr;
	if (rank == 0) {
		result = new double* [config.getSizeOfMatrix()];
		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			result[i] = new double[config.getSizeOfMatrix()];
		}
	}

	MPI_Gather(local_c[0], chunk_size * config.getSizeOfMatrix(), MPI_DOUBLE, result[0], chunk_size * config.getSizeOfMatrix(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for (long j = config.getSizeOfMatrix() - remainder; j < config.getSizeOfMatrix(); j++) {
			for (long i = 0; i < config.getSizeOfMatrix(); i++) {
				result[j][i] = matrix_a[j][i] / matrix_b[j][i];
			}
		}

		// Print the result
		cout << "MPI operation division is ended" << endl;
		double time = (double)(clock() - start_sync) / CLOCKS_PER_SEC;
		cout << "Execution time: " << time << " seconds\n" << endl;

		for (int i = 0; i < config.getSizeOfMatrix(); i++) {
			delete[] result[i];
		}
		delete[] result;
	}

	for (int i = 0; i < chunk_size; i++) {
		delete[] local_a[i];
		delete[] local_b[i];
		delete[] local_c[i];
	}
	delete[] local_a;
	delete[] local_b;
	delete[] local_c;
}