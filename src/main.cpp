#include "rotation.hpp"
#include "basicAlg.hpp"


#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"

using namespace std;

int main()
{
    int L=4, kindex, k;
    int num, numStates_m; 
    double overlap = 0, dot = 0, item, projSize;
    vector<string> states;
    complex<double> * A;
    cout.precision(17);  
    
    vector<double> groundState;
//    vector<double> groundState{
//        0,0,0,1/sqrt(6.),0,
//        0,0,1/sqrt(12.),0,1/sqrt(3.0),
//        0,0,0,0,0,
//        1/sqrt(6.), 0,  1/sqrt(12.), 0};


//    vector<double> groundState{
//        0,0,0,0,0,
//        0,0,1/sqrt(6.),0,1/sqrt(3.0),
//        0,0,0,0,0,
//        1/sqrt(3.), 0,  1/sqrt(6.), 0};


    

    ifstream myFile{"site10.txt"};
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
    
    A = new complex<double>[numStates_m]; // A is the vector for projection
    for (int i = 0; i < numStates_m; i++)
    {  A[i] = 0;}

    
    numStates_m = outputmStates(L, 0, states);

    for (k = 0; k < L; k++)
    {
        projection(A, groundState, states, k);
        projSize = normalize(A, numStates_m);
        overlap = dotProduct(A, groundState, numStates_m)/projSize/projSize;
        cout << "k = "<< k << " overlap : " << overlap << endl;
    }
    return 0;
} //  end Main




