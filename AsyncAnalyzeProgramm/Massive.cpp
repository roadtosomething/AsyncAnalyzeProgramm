#include "Massive.h"
#include "Config.h"
#include <fstream>
#include <iostream>
#include <time.h>

using namespace std;


ifstream in;
Config config;
int threads = config.getNumberOfThreads();
clock_t start_sync;
clock_t end_sync;
double seconds;

void Massive::fillingMAssive(int* array) {
	cout << "Process filling massive is started.\n";
	string fileName = config.getDataFileName();
	in.open(fileName);
	if (in.is_open()) {
		cout << "File \"" << fileName << "\" is open\n";
		cout << "Filling array...\n";
		for (long i = 0; i < config.getSizeOfArray(); i++) {
			in >> array[i];
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