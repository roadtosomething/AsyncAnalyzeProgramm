#include "Massive.h"
#include "Config.h"
#include <fstream>
#include <iostream>

using namespace std;

void Massive::fillingMAssive(int* array) {
	cout << "Process filling massive is started.\n";
	ifstream in;
	Config config;
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