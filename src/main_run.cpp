
#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"
#include "rotation.hpp"

using namespace std;



int main()
{
    int L=4, kindex, k = 0;
    int num, numStates_m; 
    double overlap = 0, dot = 0, item, projSize, groundNorm;
    vector<string> states;
    cout.precision(17);  
    
    vector<double> groundState;

    ifstream myFile{"site6_m0_open.txt"};
    int totalStates = 0;
    numStates_m = 19;
    if (!myFile.is_open()) return -1;
    myFile >> L ;
    cout << "L is : " << L << endl;
    myFile >> numStates_m ;
    cout << "numStates_m : " << numStates_m << endl;
    for (int i = 0 ; i < numStates_m;i++)
    {
        myFile >> item;
        groundState.push_back(item);

    }
    cout << "finish reading file " << endl;   
    myFile.close();
    
    groundNorm =   normalize(groundState, numStates_m) ;
    cout << "groundNorm = " << groundNorm << endl;
    for (int i =0; i < numStates_m; i++) {groundState[i] = groundState[i]/groundNorm;}



    
//    {
//    complex<double> **A = new complex<double>*[numStates_m];
//    for(int i = 0; i < numStates_m; i++)
//    {
//        A[i] = new complex<double>[numStates_m];
//    }
//    for (int i = 0; i < numStates_m; i++)
//    {
//        for (int j = 0; j < numStates_m; j++) {A[i][j] = 0;}
//    }
//    cout << "This is OK!" << endl;
    numStates_m = outputmStates(L, 0, states);
    cout << " numStates_m = " << numStates_m << endl;
    k = 0;
    for (k = 0 ; k < L; k++)
    {
        overlap = overlapKmode(groundState, k, numStates_m, states);
        cout << "k = " << k <<  " :overlap ^2: " << overlap << endl;
    }
    return 0;
} //  end Main


