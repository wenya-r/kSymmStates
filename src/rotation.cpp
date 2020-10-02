#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"

using namespace std;

void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect, vector<int> &indexTable);
int bitToNum(string bit);
string translation(string bit);

int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<int> &indexTable);
void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k, vector<int> &indexTable);
void projection(complex<double> *arr, complex<double> *groundstate , vector<string> &states, int k, vector<int> &indexTable);
int bitToIndex(string bit, vector<string> &vect, vector<int> &indexTable);

double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states, vector<int> &indexTable) ;

double overlapKmode(complex<double> *groundState, int k, int numStates_m, vector<string> states, vector<int> &indexTable) ;



/*
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
    
//    }
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

*/


double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states, vector<int> &indexTable) 
{
    double result, projSize;        

    complex<double> * A;


    A = new complex<double>[numStates_m]; // A is the vector for projection
    for (int i = 0; i < numStates_m; i++)
    {  A[i] = 0;}
//    cout << "beginning of overlpaLmode " << endl;

    projection(A, groundState, states, k, indexTable);
//    cout << "print out A " << endl;
//    for (int i = 0; i < numStates_m; i++){ cout << A[i] << endl;}

    projSize = normalize(A, numStates_m);
//    cout << "projSize = " << projSize << endl;
    result = dotProduct(A, groundState, numStates_m)/projSize/projSize;

    return result;
}

double overlapKmode(complex<double> *groundState, int k, int numStates_m, vector<string> states, vector<int> &indexTable) 
{
    double result, projSize;        

    complex<double> * A;


    A = new complex<double>[numStates_m]; // A is the vector for projection
    for (int i = 0; i < numStates_m; i++)
    {  A[i] = 0;}
//    cout << "beginning of overlpaLmode " << endl;

    projection(A, groundState, states, k, indexTable);


    projSize = normalize(A, numStates_m);
    cout << endl;
    result = dotProduct(A, groundState, numStates_m)/projSize/projSize;

    return result;
}




void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k, vector<int> &indexTable)
{
    int L, size = states.size();
    string state;
    const double pi = acos(-1);
    double  km;
    complex<double> phase, phaser;
    L = states[0].length() ; 
    km = k*2*pi/L;
    phase = complex<double>(0, km);
    cout << " states size: " << size << endl;
    cout << "beginning of projection vec<double> " << endl;
    for (int i = 0; i< size; i++)
    {
        arr[i] = groundstate[i];
        state = states[i];
//        cout << "i = " << i << arr[i] << endl;
        phaser= phase;
        for (int j = 1; j < L ; j++)
        { 
            state = translation(state);            
            arr[i] = arr[i] + exp(phaser) * groundstate[bitToIndex(state,states, indexTable)] ;
            phaser = phaser + phase;
//            cout << arr[i] << endl;
        }
//        cout << arr[i] << endl; // Problem!!
    }    
}


void projection(complex<double> *arr, complex<double> *groundstate , vector<string> &states, int k, vector<int> &indexTable)
{
    int L, size = states.size();
    string state;
    const double pi = acos(-1);
    double  km;
    complex<double> phase, phaser;
    L = states[0].length() ; 
    km = k*2*pi/L;
    phase = complex<double>(0, km);
//    cout << "size : " << size << endl;
//    cout << "beginning of projection " << endl;
    for (int i = 0; i< size; i++)
    {
        arr[i] = groundstate[i];
        state = states[i];
        phaser= phase;

//        cout << "i = " << i << endl;
//        cout << "phase = " << phaser << endl;
//        cout << arr[i] << endl;
        for (int j = 1; j < L ; j++)
        { 
            state = translation(state);            
            arr[i] = arr[i] + exp(phaser) * groundstate[bitToIndex(state,states, indexTable)] ;
            phaser = phaser + phase;
//            cout << arr[i] << endl;
        }
//        cout << arr[i] << endl;
    }    
}

void outputAllStates(int L){
    int numOfStates = pow(3,L);
    int stateNum;
    string state;
    cout << "num of Allstates: " << numOfStates << endl;
    for(int i = numOfStates-1;i>-1;i--)
    {   
//        cout << i;
        stateNum = i;
        state = to_string(stateNum%3);
        for (int j = 0; j < 3; j++)
        {
            stateNum = stateNum/3;
            state = to_string(stateNum%3) + state;
        }
//        cout << state << endl;
    }
}

int outputmStates(int L, int m, vector<string> &vect, vector<int> &indexTable){
// m is Sz total. sum(string)= L+m 

    int numOfStates = pow(3,L);
    int stateNum, bit, sumBit;
    int total = 0;
    string state;
    cout << "num of m states: " << numOfStates << endl;
    for(int i = numOfStates-1;i>-1;i--)
    {   
        stateNum = i;
        
        bit = stateNum%3;
        state = to_string(bit);
        sumBit = bit;
        for (int j = 0; j < L-1; j++)
        {
            stateNum = stateNum/3;
            bit =  stateNum%3;
            state = to_string(bit) + state;
            sumBit = sumBit + bit;
        }
        if (sumBit == L+m)
        {
            total++;
            vect.push_back(state);
            indexTable.push_back(i);
//            cout << state << endl;
        }
    }
    return total;

}

int bitToNum(string bit)
{
    int L, num;
    num = 0;
//    cout << "do you see bit?" << bit<< endl;
    L = bit.length();
    num = stoi(bit, nullptr, 3);
    return num;
}



string translation(string bit)
{
    string transBit ;
    transBit = bit[1];
    for (int i = 2; i < bit.length(); i++)
    {
        transBit = transBit+bit[i];
    }
    transBit = transBit + bit[0];
    return transBit;

}

int bitToIndex(string bit, vector<string> &vect, vector<int> &indexTable)
{
    int numBit, index= 0, upper,lower, mid=0;
    bool done= false;
    numBit = bitToNum(bit);
//    cout << "numBit = " << numBit<< endl;
//    while(!done)
//    {
//        if (vect[index] == bit)
//        {
//            done = true;
//        }
//        else{index++;}
//    }
//
//    for (int j : indexTable){cout << j << "j" << endl;}
    upper = 0; 
    lower = vect.size()-1;
//    cout << "inside bitToIndex " << numBit<< endl;
    while(lower > upper || lower == upper)
    {
        mid = (upper + lower)/2;
//        cout << "upper " << upper << endl;
//        cout << "lower " << lower << endl;
//        cout << " mid  " << mid   << endl;
//        cout << "index[]" << indexTable[mid] << endl;
        if (indexTable[mid]> numBit){  upper = mid+1 ; }
        else if (indexTable[mid]< numBit) { lower = mid-1;}
        else{return mid;}
    }
    return mid;
}




int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<int> &indexTable)
{
    string bit;
    int i, r, kindex = 0;
    double km, norm, product ;
    bool iskState;
    const double pi = acos(-1);
    complex<double> * cofArr;
    const complex<double> imi(0,1);
    complex<double> phase, phaser;
  
    //declare complex array and set to zero
    cofArr = new complex<double>[col];
    for (int i = 0; i < col; i++) {cofArr[i]=0;}

    km = k*2*pi/L;
    phase = complex<double>(0, km);
    for(i = 0; i < row; i++)
    {
        bit = vect[i];
        for (int i = 0; i < col; i++) {cofArr[i]=0;}
        cofArr[bitToIndex(bit, vect, indexTable)] = 1;        
        phaser = phase;
        for(r = 1; r < L; r++)
        {
            bit = translation(bit);
            cofArr[bitToIndex(bit, vect, indexTable)] =  cofArr[bitToIndex(bit, vect, indexTable)]+exp(phaser);   
            phaser = phaser + phase;
        }
        norm =  normalize(cofArr, col) ;

    // normalize the vector
        if (norm > 0.00000001)
        {
            for (r = 0; r < col; r++)
            {
                cofArr[r] = cofArr[r]/norm;
            }
            // check orthoganal with all previous
            iskState = true;
            for (int j = 0; j <kindex; j++)
            {
                if (dotProduct(cofArr, arr[j], col) > 0.5){iskState = false;}
            }
            if(iskState)
            {
                for (int j=0; j < col; j++)
                { arr[kindex][j] = cofArr[j];} //  to do
                kindex++;
            }
        }
    }       
    return kindex;
}
