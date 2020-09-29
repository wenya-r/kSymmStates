#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"

using namespace std;

void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect);
int bitToNum(string bit);
string translation(string bit);

int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect);
void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k);
int bitToIndex(string bit, vector<string> &vect);

double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states) ;

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


double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states) 
{
    double result, projSize;        

    complex<double> * A;


    A = new complex<double>[numStates_m]; // A is the vector for projection
    for (int i = 0; i < numStates_m; i++)
    {  A[i] = 0;}


    projection(A, groundState, states, k);


    projSize = normalize(A, numStates_m);
    result = dotProduct(A, groundState, numStates_m)/projSize/projSize;

    return result;
}

void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k)
{
    int L, size = states.size();
    string state;
    const double pi = acos(-1);
    double  km;
    complex<double> phase, phaser;
    L = states[0].length() ; 
    km = k*2*pi/L;
    phase = complex<double>(0, km);
    for (int i = 0; i< size; i++)
    {
        arr[i] = groundstate[i];
        state = states[i];
        phaser= phase;
        for (int j = 1; j < L ; j++)
        { 
            state = translation(state);            
            arr[i] = arr[i] + exp(phaser) * groundstate[bitToIndex(state,states)] ;
            phaser = phaser + phase;
        }
        
    }    



}



void outputAllStates(int L){
    int numOfStates = pow(3,L);
    int stateNum;
    string state;
    cout << "num of states: " << numOfStates << endl;
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

int outputmStates(int L, int m, vector<string> &vect){
// m is Sz total. sum(string)= L+m 

    int numOfStates = pow(3,L);
    int stateNum, bit, sumBit;
    int total = 0;
    string state;
    cout << "num of states: " << numOfStates << endl;
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

int bitToIndex(string bit, vector<string> &vect)
{
    int index= 0;
    bool done= false;
    while(!done)
    {
        if (vect[index] == bit)
        {
            done = true;
        }
        else{index++;}
    }

    return index;
}




int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect)
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
        cofArr[bitToIndex(bit, vect)] = 1;        
        phaser = phase;
        for(r = 1; r < L; r++)
        {
            bit = translation(bit);
            cofArr[bitToIndex(bit, vect)] =  cofArr[bitToIndex(bit, vect)]+exp(phaser);   
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



//
//double dotProduct(complex<double> * vec1, complex<double> * vec2, int col)
//{
//    complex<double> product = 0;
//    for (int i = 0; i < col; i++)
//    {        product = product + vec1[i]*conj(vec2[i]); }
//    return norm(product);
//}
//
//
//double dotProduct(complex<double> * vec1, vector<double> vec2, int col)
//{
//    complex<double> product = 0;
//    for (int i = 0; i < col; i++)
//    {        product = product + vec1[i].real()*(vec2[i]) + complex<double>(0,1)*vec1[i].imag()*(vec2[i]); }
//    return norm(product); //norm : squared magnitude - The norm of (3,4) is 25
//}


