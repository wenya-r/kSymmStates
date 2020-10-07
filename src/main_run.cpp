
#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"
#include "rotation.hpp"
#include "FileIO.hpp"

using namespace std;



int main()
{
    int L=22, k = 0;
    long long int ni, sz, num, numStates_m, numStates; 
    double overlap = 0, dot = 0, item, projSize, groundNorm;
    vector<string> states;
    cout.precision(17);  
    FILE *fp;
//    vector<double> groundState;
    vector<unsigned long int> indexTable;
    complex<double> arr;
    complex<double> * groundState;


    numStates_m = outputmStates(L, 0, indexTable);
    fp = fopen("site22zvo_eigenvec_0_rank_0.dat", "rb");
        if(fp == NULL){
        printf("file error");
        exit(2);
    }
    cout << "sizeof(ni) " << sizeof(numStates) << endl;
    fread(&ni, 4, 1, fp);
    fread(&numStates, 10, 1, fp);
    printf("number of iterations = %d\n",ni );
    printf("number of m_states from file = %d\n",numStates );
    
//    fp.close();

    groundState = new complex<double>[numStates_m]; // A is the vector for projection
    fread(&arr, sizeof(arr), 1, fp);
    for (long int i = 0; i < numStates_m; i++)
    {  
        fread(&arr, sizeof(arr), 1, fp);
//        cout << sqrt(norm(arr)) << endl;
        groundState[i] = arr;
        while (i < 50)
        {cout << arr << endl;}
    }
    cout << "L is : " << L << endl;
//    myFile >> numStates_m ;
//    cout << "numStates_m : " << numStates_m << endl;
//    for (int i = 0 ; i < numStates_m;i++)
//    {
//        myFile >> item;
//        groundState.push_back(item);
//
//    }
    cout << "finish reading file " << endl;   
//    myFile.close();
    
    groundNorm =   dotProduct(groundState, groundState, numStates_m) ;
    cout << "groundNorm = " << groundNorm << endl;
    for (int i =0; i < numStates_m; i++) {groundState[i] = groundState[i]/groundNorm;}


    // save indexTable
    saveStates(indexTable, "indexTable20site.dat");
    
//    readStates(indexTable, "indexTable22site.dat");

    cout << " numStates_m = " << numStates_m << endl;
    k = 9;
    for (k = 0; k < L; k++)
    {
        overlap = overlapKmode(L, groundState, k, numStates_m, indexTable);
        cout << "k = " << k <<  " :overlap ^2: " << overlap << endl;
    }
    return 0;
} //  end Main


