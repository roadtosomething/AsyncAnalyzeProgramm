#pragma once
#include "Config.h"
#include <time.h>
#include <iostream>

class Massive
{
public:
	void fillingMAssive(int*);
	unsigned long sumOfArray(int*);
	unsigned long sumOfArrayAsync(int*);
	unsigned long sumOfArrayAsyncMPI(int*);
	void quicksort(int*, long, long);
	void quicksortAsync(int*, long, long);
private:
	long partition(int*, long, long);
	long partitionAsync(int*, long, long);
	Config config;
	int threads = config.getNumberOfThreads();
	clock_t start_sync;
	clock_t end_sync;
	double seconds;
};