#pragma once 

#include <vector>
#include <string>
#include <complex>


void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect);
int bitToNum(string bit);
string translation(string bit);

int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect);
void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k);
int bitToIndex(string bit, vector<string> &vect);

double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states) ;


