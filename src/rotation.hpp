#pragma once 

#include <vector>
#include <string>
#include <complex>


void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect, vector<unsigned long int> &indexTable);
unsigned long int bitToNum(string bit);
string translation(string bit);

int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<unsigned long int> &indexTable);
void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k, vector<unsigned long int> &indexTable);
int bitToIndex(string bit, vector<string> &vect, vector<unsigned long int> &indexTable);

double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states, vector<unsigned long int> &indexTable) ;

double overlapKmode(complex<double> *groundState, int k, int numStates_m, vector<string> states, vector<unsigned long int> &indexTable) ;

