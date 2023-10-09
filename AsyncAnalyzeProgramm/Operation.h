#pragma once
#include "mpi.h"
class Operation
{
public:
	void addition(int*, int*);
	void addition(int**, int**);
	void subtraction(int*, int*);
	void subtraction(int**, int**);
	void multiplication(int*, int*);
	void multiplication(int**, int**);
	void division(int*, int*);
	void division(int**, int**);
	//OMP
	void additionAsync(int*, int*);
	void additionAsync(int**, int**);
	void subtractionAsync(int*, int*);
	void subtractionAsync(int**, int**);
	void multiplicationAsync(int*, int*);
	void multiplicationAsync(int**, int**);
	void divisionAsync(int*, int*);
	void divisionAsync(int**, int**);
	//MPI
	void additionMPI(int*, int*);
	void additionMPI(int**, int**);
	void subtractionMPI(int*, int*);
	void subtractionMPI(int**, int**);
	void multiplicationMPI(int*, int*);
	void multiplicationMPI(int**, int**);
	void divisionMPI(int*, int*);
	void divisionMPI(int**, int**);
private:
	MPI_Status status;
	int rank, size;
};

