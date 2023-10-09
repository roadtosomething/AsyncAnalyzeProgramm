#include "TestingData.h"
#include <fstream>


void TestingData::crateNewDataFile(string fileName){
    ofstream out;
    long i, N;
    N = 5000000;
    out.open(fileName);
    for (i = 0; i < N; i++)
    {
        if (out.is_open()) {
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << " ";
            out << rand() << std::endl;
        }
        if ((i + 1) % (N / 100) == 0) { std::cout << "Task 1 completed: " << (i + 1) / (N / 100) << "/100%\n"; }
    }
    out.close();
}