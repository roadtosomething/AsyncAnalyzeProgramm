#pragma once
#include <iostream>

using namespace std;

class Config
{
private:
	int number_of_threads = 10;
	long size_of_array = 1000000;
	long size_of_matrix = 1000;
	string data_file_name = "data.txt";
public:
	int getNumberOfThreads();
	long getSizeOfArray();
	long getSizeOfMatrix();
	string getDataFileName();
};

