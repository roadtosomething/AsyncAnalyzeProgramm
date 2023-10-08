#pragma once
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
	void additionAsync(int*, int*);
	void additionAsync(int**, int**);
	void subtractionAsync(int*, int*);
	void subtractionAsync(int**, int**);
	void multiplicationAsync(int*, int*);
	void multiplicationAsync(int**, int**);
	void divisionAsync(int*, int*);
	void divisionAsync(int**, int**);
};

