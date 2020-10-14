
#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"
#include "rotation.hpp"
#include "FileIO.hpp"
#include <chrono>
#include <set>


using namespace std;

set<int> indexT;

int main(int argc, char **argv) // L , filename, SzTotal, k
{
    int L, k, SzTotal = 0;
    long long int ni, sz, num, numStates_m, numStates; 
    double overlap = 0, dot = 0, item, projSize, groundNorm;
    vector<string> states;
    cout.precision(17);  
    FILE *fp;
    vector<unsigned long int> indexTable;
    complex<double> arr;
    complex<double> * groundState;
    string fileName;
    string indexName;

//    cout << argv[4] << endl;
    L = stoi(argv[1]);
    cout << "L = " << L << endl;
    fileName = argv[2];
    cout << "L = " << fileName << endl;
    SzTotal = stoi(argv[3]);
    cout << "m = " << SzTotal << endl;
    k = stoi(argv[4]);
    fp = fopen(fileName.c_str(), "rb");
        if(fp == NULL){
        printf("file error");
        exit(2);
    }
    cout << "sizeof(ni) " << sizeof(numStates) << endl;
    cout << "sizeof(arr) " << sizeof(arr) << endl;
    fread(&ni, 4, 1, fp);
    fread(&numStates, 8, 1, fp);
    printf("number of iterations = %d\n",ni );
    printf("number of m_states from file = %d\n",numStates );
    
//    fp.close();

    groundState = new complex<double>[numStates]; // A is the vector for projection
    fread(&arr, 16, 1, fp);
    for (int i = 0; i < numStates; i++)
    {  
        fread(&arr, sizeof(arr), 1, fp);
        groundState[i] = arr;
    }
    cout << "L is : " << L << endl;
    cout << "finish reading file " << endl;   
    auto start = chrono::high_resolution_clock::now();
    numStates_m = outputmStates(L, SzTotal, indexTable);
//    numStates_m = outputmStates(L, SzTotal, indexT);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "time spent gen. m in seconds :" << duration.count() << endl;


    groundNorm =   dotProduct(groundState, groundState, numStates) ;
    cout << "groundNorm = " << groundNorm << endl;
    for (int i =0; i < numStates; i++) {groundState[i] = groundState[i]/groundNorm;}

//  conver indexT to indexTable
//    copy(indexT.begin(), indexT.end(), back_inserter(indexTable));
//    reverse(indexTable.begin(), indexTable.end());

//  end converting




    indexName = "indexTableL"+to_string(L)+"m"+to_string(SzTotal)+".dat";
    // save indexTable
    saveStates(indexTable, indexName);
    
//    readStates(indexTable, indexName);

    cout << " numStates from indexfile= " << indexTable.size() << endl;

    for (k = k; k < L; k++)
    {
        overlap = overlapKmode(L, groundState, k, numStates, indexTable);
        cout << "k = " << k <<  " :overlap ^2: " << overlap << endl;
    }
    return 0;
} //  end Main


