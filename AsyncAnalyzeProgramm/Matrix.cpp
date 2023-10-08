#include "Matrix.h"

#include "Config.h"
#include <fstream>
#include <iostream>

using namespace std;


void Matrix::fillingMatrix (int** matrix) {
	cout << "Process filling matrix is started.\n";
	ifstream in;
	Config config;
	string fileName = config.getDataFileName();
	in.open(fileName);
	if (in.is_open()) {
		cout << "File \"" << fileName << "\" is open\n";
		cout << "Filling array...\n";
		for (long i = 0; i < config.getSizeOfMatrix(); i++) {
			for (long j = 0; j < config.getSizeOfMatrix(); j++) {
				in >> matrix[i][j];
			}
		}
		cout << "Array is filling.\n";
	}
	else {
		cout << "File \"" << fileName << "\" is not open. Check if this file exists.\n";
	}
	cout << "Process filling matrix is ended.\n";
}