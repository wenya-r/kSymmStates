#pragma once 

#include <vector>
#include <string>
#include <complex>


void outputAllStates(int L);
int outputmStates(int L, int m, vector<unsigned long int> &indexTable);
unsigned long int bitToNum(string bit);
string translation(string bit);

int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<unsigned long int> &indexTable);
void projection(int L, complex<double> *arr, vector<double> groundstate , int k, vector<unsigned long int> &indexTable);
void projection(int L, complex<double> *arr, complex<double> *groundstate, int k, vector<unsigned long int> &indexTable);
int bitToIndex(string bit, vector<unsigned long int> &indexTable);

double overlapKmode(int L, vector<double> groundState, int k, int numStates_m, vector<unsigned long int> &indexTable) ;

double overlapKmode(int L, complex<double> *groundState, int k, int numStates_m, vector<unsigned long int> &indexTable) ;

string numToBit(int L, unsigned long int num);
