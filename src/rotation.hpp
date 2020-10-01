#pragma once 

#include <vector>
#include <string>
#include <complex>


void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect, vector<int> &indexTable);
int bitToNum(string bit);
string translation(string bit);

int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<int> &indexTable);
void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k, vector<int> &indexTable);
int bitToIndex(string bit, vector<string> &vect, vector<int> &indexTable);

double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states, vector<int> &indexTable) ;

double overlapKmode(complex<double> *groundState, int k, int numStates_m, vector<string> states, vector<int> &indexTable) ;

