#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>
#include "basicAlg.hpp"
#include "rotation.hpp"

using namespace std;

//void outputAllStates(int L);
//int outputmStates(int L, int m, vector<string> &vect, vector<unsigned long int> &indexTable);
//unsigned long int bitToNum(string bit);
//
//string translation(string bit);
//
//int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<unsigned long int> &indexTable);
//void projection(complex<double> *arr, vector<double> groundstate , vector<string> &states, int k, vector<unsigned long int> &indexTable);
//void projection(complex<double> *arr, complex<double> *groundstate , vector<string> &states, int k, vector<unsigned long int> &indexTable);
//int bitToIndex(string bit, vector<string> &vect, vector<unsigned long int> &indexTable);

//double overlapKmode(vector<double> groundState, int k, int numStates_m, vector<string> states, vector<unsigned long int> &indexTable) ;
//
//double overlapKmode(complex<double> *groundState, int k, int numStates_m, vector<string> states, vector<unsigned long int> &indexTable) ;





double overlapKmode(int L, vector<double> groundState, int k, int numStates_m, vector<unsigned long int> &indexTable) 
{
    double result, projSize;        

    complex<double> * A;


    A = new complex<double>[numStates_m]; // A is the vector for projection
    for (int i = 0; i < numStates_m; i++)
    {  A[i] = 0;}
    cout << "beginning of overlpaLmode " << endl;

    projection(L, A, groundState, k, indexTable);
//    cout << "print out A " << endl;
//    for (int i = 0; i < numStates_m; i++){ cout << A[i] << endl;}

    projSize = normalize(A, numStates_m);
//    cout << "projSize = " << projSize << endl;
    result = dotProduct(A, groundState, numStates_m)/projSize/projSize;

    return result;
}

double overlapKmode(int L, complex<double> *groundState, int k, int numStates_m, vector<unsigned long int> &indexTable) 
{
    double result, projSize;        

    complex<double> * A;


    A = new complex<double>[numStates_m]; // A is the vector for projection
    for (int i = 0; i < numStates_m; i++)
    {  A[i] = 0;}
    cout << "beginning of overlpaLmode " << endl;

    projection(L, A, groundState, k, indexTable);


    projSize = normalize(A, numStates_m);
    cout << endl;
    result = dotProduct(A, groundState, numStates_m)/projSize/projSize;

    return result;
}




void projection(int L, complex<double> *arr, vector<double> groundstate , int k, vector<unsigned long int> &indexTable)
{
    int size = groundstate.size();
    string state;
    const double pi = acos(-1);
    double  km;
    complex<double> phase, phaser;
    L = 4; //states[0].length() ; 
    km = k*2*pi/L;
    phase = complex<double>(0, km);
//    cout << " states size: " << size << endl;
    cout << "beginning of projection vec<double> " << endl;
    for (int i = 0; i< size; i++)
    {
//        cout << "begin loop" << arr[i] << groundstate[i] << endl;
        arr[i] = groundstate[i];


       state = numToBit(L, indexTable[i]);

//        cout << "i = " << i << state<< endl;
        phaser= phase;
        for (int j = 1; j < L ; j++)
        { 
            state = translation(state);            
            arr[i] = arr[i] + exp(phaser) * groundstate[bitToIndex(state, indexTable)] ;
            phaser = phaser + phase;
//            cout << arr[i] << endl;
        }
//        cout << arr[i] << endl; 
    }    
}


void projection(int L, complex<double> *arr, complex<double> *groundstate, int k, vector<unsigned long int> &indexTable)
{
    int size = indexTable.size();
    string state;
    const double pi = acos(-1);
    double  km;
    complex<double> phase, phaser;
    km = k*2*pi/L;
    phase = complex<double>(0, km);
    cout << "size : " << size << endl;
    cout << "beginning of projection " << endl;
    for (int i = 0; i< size; i++)
    {
        arr[i] = groundstate[i];
//        state = states[i];
        state = numToBit(L, indexTable[i]);
        phaser= phase;

        cout << "i = " << state << endl;
//        cout << "phase = " << phaser << endl;
//        cout << arr[i] << endl;
        for (int j = 1; j < L ; j++)
        { 
            state = translation(state);            
            arr[i] = arr[i] + exp(phaser) * groundstate[bitToIndex(state, indexTable)] ;


            phaser = phaser + phase;
//            cout << arr[i] << endl;
        }
//        cout << arr[i] << endl;
    }    
}

void outputAllStates(int L){
    unsigned long int numOfStates = pow(3,L);
    unsigned long int stateNum;
    string state;
    cout << "num of Allstates: " << numOfStates << endl;
//    for(unsigned long int i = numOfStates-1;i>-1;i--)
//    {   
//        cout << i;
//        stateNum = i;
//        state = to_string(stateNum%3);
//        for (int j = 0; j < 3; j++)
//        {
//            stateNum = stateNum/3;
//            state = to_string(stateNum%3) + state;
//        }
//        cout << state << endl;
//    }
}

int outputmStates(int L, int m, vector<unsigned long int> &indexTable){
// m is Sz total. sum(string)= L+m 

    long int numOfStates = pow(3,L);
    long int stateNum, bit, sumBit;
    long int i;

    int total = 0;
    string state;

    printf(" L = %d \n " , L );
    cout <<  "num of m states: " ;
    printf("%ld\n", numOfStates);
    i = numOfStates - 1 ;
//    cout << "max_size of vector indexTable = " << indexTable.max_size() << endl;
//    cout << " i is : " << (i > -1) << endl ;
//    printf("%ld\n", i);
    while(i > -1)
    {   

//    cout << " i is : " ;
//    printf("%ld\n", i);

//        cout << " inside i " << endl;
        stateNum = i;
        
        bit = stateNum%3;
        state = to_string(bit);
        sumBit = bit;
        if (i > 3486784380)           
        {
            cout << "stateNum" << stateNum << endl;
        }


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
//            vect.push_back(state);
            indexTable.push_back(i);
//            cout << state << endl;
        }
        i--;
    }
    

    return total;

}

unsigned long int bitToNum(string bit)
{
    unsigned long int L, num;
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

int bitToIndex(string bit, vector<unsigned long int> &indexTable)
{
    int index= 0, upper,lower, mid=0;
    unsigned long int numBit;
    bool done= false;
    numBit = bitToNum(bit);
    upper = 0; 
    lower = indexTable.size()-1;
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




int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect, vector<unsigned long int> &indexTable)
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
        cofArr[bitToIndex(bit, indexTable)] = 1;        
        phaser = phase;
        for(r = 1; r < L; r++)
        {
            bit = translation(bit);
            cofArr[bitToIndex(bit, indexTable)] =  cofArr[bitToIndex(bit, indexTable)]+exp(phaser);   
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



string numToBit(int L, unsigned long int num)
{
    string item;
    item = to_string(num%3);
    for (int i = 1; i < L ; i++)
    {
        num = num / 3;
        item = to_string(num%3) + item; 

    }
    return item;


}












