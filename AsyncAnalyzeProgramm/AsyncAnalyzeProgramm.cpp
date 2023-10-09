
#include "Config.h"
#include "Operation.h"
#include "Massive.h"
#include "Matrix.h"
#include "TestingData.h"

using namespace std;



int main(int argc, char** argv)
{
	//initialize parameters
	Config config;
	long arraySize = config.getSizeOfArray();
	long matrixSize = config.getSizeOfMatrix();
	TestingData testingData;
	Operation operation;
	Massive massive;
	Matrix matrix;
	clock_t start_sync;
	clock_t end_sync;
	double seconds;
	int* array_A = new int[arraySize];

	int* array_b = new int[arraySize];

	int** matrix_A = new int*[matrixSize];
	for (long i = 0; i < matrixSize; i++) {
		matrix_A[i] = new int[matrixSize];
	}

	int** matrix_B = new int* [matrixSize];
	for (long i = 0; i < matrixSize; i++) {
		matrix_B[i] = new int[matrixSize];
	}

	//Start programm
	cout << "Start Analyze programm" << endl;
	//Create file with data
	cout << "Need a new data file for filling matrix and arrays?\n 1. Yes\n 2. No" << endl;
	short answer=2;
	cin >> answer;
	if (answer == 1) {
		cout << "Start creating a new file with name \"" << config.getDataFileName() << "\"" << endl;
		testingData.crateNewDataFile(config.getDataFileName());
		cout << "File is created" << endl;
	}
	massive.fillingMAssive(array_A);
	massive.fillingMAssive(array_b);
	matrix.fillingMatrix(matrix_A);
	matrix.fillingMatrix(matrix_B);
	//Task 1
	cout << "Task1\n" << endl;
	massive.sumOfArray(array_A);
	massive.sumOfArrayAsync(array_A);
	cout << "|--------------------|"<<endl;
	//Task 2
	cout << "Task 2\n" << endl;
	cout << "Process of sync quicksort is started..." << endl;
	start_sync = clock();
	massive.quicksort(array_A, 0, arraySize-1);
	end_sync = clock();
	cout << "Process of quicksort is ended." << endl;
	seconds = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << seconds << " seconds\n" << endl;
	cout << "Process of async quicksort is started..." << endl;
	start_sync = clock();
	massive.quicksortAsync(array_b, 0, arraySize-1);
	end_sync = clock();
	cout << "Process of async quicksort is ended." << endl;
	seconds = (double) (end_sync - start_sync) / CLOCKS_PER_SEC;
	cout << "Execution time: " << seconds << " seconds\n" << endl;
	cout << "|--------------------|" << endl;
	// Task 3
	cout << "Task 3\n" << endl;
	cout << "Synchronize" << endl;
	operation.addition(array_A, array_b);
	operation.subtraction(array_A, array_b);
	operation.multiplication(array_A, array_b);
	operation.division(array_A, array_b);
	cout << "|--------------------|" << endl;
	cout << "OMP" << endl;
	operation.additionAsync(array_A, array_b);
	operation.subtractionAsync(array_A, array_b);
	operation.multiplicationAsync(array_A, array_b);
	operation.divisionAsync(array_A, array_b);
	cout << "|--------------------|" << endl;
	//Task 4
	cout << "Task 4\n" << endl;
	cout << "Synchronize" << endl;
	operation.addition(matrix_A, matrix_B);
	operation.subtraction(matrix_A, matrix_B);
	operation.multiplication(matrix_A, matrix_B);
	operation.division(matrix_A, matrix_B);
	cout << "|--------------------|" << endl;
	cout << "OMP" << endl;
	operation.additionAsync(matrix_A,matrix_B);
	operation.subtractionAsync(matrix_A,matrix_B);
	operation.multiplicationAsync(matrix_A,matrix_B);
	operation.divisionAsync(matrix_A,matrix_B);
	cout << "|--------------------|" << endl;
	cout << "MPI" << endl;
	cout << "Task 3" << endl;
	MPI_Init(&argc, &argv);
	operation.additionMPI(array_A, array_b);
	operation.subtractionMPI(array_A, array_b);
	operation.multiplicationMPI(array_A, array_b);
	operation.divisionMPI(array_A, array_b);
	cout << "|--------------------|" << endl;
	cout << "Task 4" << endl;
	operation.additionMPI(matrix_A, matrix_B);
	operation.subtractionMPI(matrix_A, matrix_B);
	operation.multiplicationMPI(matrix_A, matrix_B);
	operation.divisionMPI(matrix_A, matrix_B);
	cout << "|--------------------|" << endl;
	MPI_Finalize();
	//Deleted object
	delete[] array_A;
	delete[] array_b;
	for (long i = 0; i < matrixSize; i++) {
		delete[] matrix_A[i];
	}
	delete[] matrix_A;
	for (long i = 0; i < matrixSize; i++) {
		delete[] matrix_B[i];
	}
	delete[] matrix_B;
}