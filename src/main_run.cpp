
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
    int L=20, k = 0, ni, sz;
    int num, numStates_m; 
    double overlap = 0, dot = 0, item, projSize, groundNorm;
    vector<string> states;
    cout.precision(17);  
    FILE *fp;
//    vector<double> groundState;
    vector<unsigned long int> indexTable;
    complex<double> arr;
    complex<double> * groundState;


    fp = fopen("site20zvo_eigenvec_0_rank_0.dat", "rb");
        if(fp == NULL){
        printf("file error");
        exit(2);
    }

    fread(&ni, 4, 1, fp);
    fread(&numStates_m, 8, 1, fp);
    printf("number of iterations = %d\n",ni );
    printf("number of m_states from file = %d\n",numStates_m );
    
//    fp.close();

    groundState = new complex<double>[numStates_m]; // A is the vector for projection
    fread(&arr, sizeof(arr), 1, fp);
    for (int i = 0; i < numStates_m; i++)
    {  
        fread(&arr, sizeof(arr), 1, fp);
//        cout << sqrt(norm(arr)) << endl;
        groundState[i] = arr;
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
//    groundNorm =   normalize(groundState, numStates_m) ;
    cout << "groundNorm = " << groundNorm << endl;
    for (int i =0; i < numStates_m; i++) {groundState[i] = groundState[i]/groundNorm;}


    numStates_m = outputmStates(L, 0, indexTable);
    // save indexTable
    //
    ofstream wf("indexTable.dat", ios::out | ios::binary);
    if(!wf)
    {
        cout << "Cannot open file!" << endl;
        return 1;
    }
    for (int i = 0; i < numStates_m; i++)
        wf.write((char *) &indexTable[i], sizeof(indexTable[0]));
    wf.close();
    if (!wf.good()) 
    {
        cout << "Error occured at writing time!" << endl;
        return 1;
    }


    



    cout << " numStates_m = " << numStates_m << endl;
    k = 0;
    for (k = 0 ; k < L; k++)
    {
        overlap = overlapKmode(L, groundState, k, numStates_m, indexTable);
        cout << "k = " << k <<  " :overlap ^2: " << overlap << endl;
    }
    return 0;
} //  end Main


