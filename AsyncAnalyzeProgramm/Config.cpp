#include "Config.h"

string Config::getDataFileName() {
	return data_file_name;
}

int Config::getNumberOfThreads() {
	return number_of_threads;
}

long Config::getSizeOfArray() {
	return size_of_array;
}

long Config::getSizeOfMatrix() {
	return size_of_matrix;
}