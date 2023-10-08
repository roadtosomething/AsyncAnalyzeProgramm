#include <iostream>
#include "Config.h"
#include "Operation.h"
#include "Massive.h"
#include "Matrix.h"
#include "TestingData.h"

using namespace std;
Config config;
long arraySize = config.getSizeOfArray();
long matrixSize = config.getSizeOfMatrix();
TestingData testingData;

int main()
{
	int* array_A = new int[arraySize];
	int* array_b = new int[arraySize];
	int** matrix_A = new int*[matrixSize];
	for (long i = 0; i < matrixSize; i++) {
		matrix_A[i] = new int[matrixSize];
	}

	cout << "Start Analyze programm" << endl;
	cout << "Need a new data file for filling matrix and arrays?\n 1. Yes\n 2. No" << endl;
	short answer=2;
	cin >> answer;
	if (answer == 1) {
		cout << "Start creating a new file with name \"" << config.getDataFileName() << "\"" << endl;
		testingData.crateNewDataFile(config.getDataFileName());
		cout << "File is created" << endl;
	}
	cout << "Task1\n" << endl;

}