#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>





using namespace std;

void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect);
int bitToNum(string bit);
string translation(string bit);

double dotProduct(complex<double> * vec1, complex<double> * vec2, int col);
double dotProduct(complex<double> * vec1, vector<double> vec2, int col);
int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect);


int main()
{
    int L, kindex;
    int num, numStates_m; 
    double overlap = 0, dot = 0, item;
    vector<string> states;
    cout.precision(17);  
    
    vector<double> groundState;

    ifstream myFile{"site6_m0_open.txt"};
    int totalStates = 0;

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

    complex<double> **A = new complex<double>*[numStates_m];
    for(int i = 0; i < numStates_m; i++)
    {
        A[i] = new complex<double>[numStates_m];
    }
    for (int i = 0; i < numStates_m; i++)
    {
        for (int j = 0; j < numStates_m; j++) {A[i][j] = 0;}
    }

    numStates_m = outputmStates(L, 0, states);
//    for (int i = 0; i < numStates_m; i++)
//    {    cout << groundState[i] << endl;}
    kindex = coefficientskMomentum(A, numStates_m, numStates_m, L, 5, states);
    cout << "kinde: " << kindex << endl;
    for (int i = 0; i < kindex; i++)
    {
        dot = dotProduct(A[i], groundState, numStates_m);

        cout<< "output each overlap:" << dot << endl;
        overlap = overlap + dot;
    }
    cout << "overlap : " << overlap << endl;
//    cout << "overlap sqrt: " << sqrt(overlap) << endl;

    return 0;
} //  end Main




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



double normalize(complex<double> * vec, int col)
{
    double norm = 0;
    for (int i = 0; i < col; i++)
    {
        norm = norm + vec[i].real() * vec[i].real() + vec[i].imag() * vec[i].imag(); 
    }

    return sqrt(norm);

}


double dotProduct(complex<double> * vec1, complex<double> * vec2, int col)
{
    complex<double> product = 0;
    for (int i = 0; i < col; i++)
    {        product = product + vec1[i]*conj(vec2[i]); }
    return norm(product);
}


double dotProduct(complex<double> * vec1, vector<double> vec2, int col)
{
    complex<double> product = 0;
    for (int i = 0; i < col; i++)
    {        product = product + vec1[i].real()*(vec2[i]) + complex<double>(0,1)*vec1[i].imag()*(vec2[i]); }
    return norm(product); //norm : squared magnitude - The norm of (3,4) is 25
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
    cout << " the number of cofAff = " << sizeof(cofArr) << endl;
//
    km = k*2*pi/L;
    
    phase = complex<double>(0, km);
    for(i = 0; i < row; i++)
    {
//        cout << i << endl;
        bit = vect[i];
        for (int i = 0; i < col; i++) {cofArr[i]=0;}
        //construct the first state
        // should be a function find out the coeff for |a> and k 
        
        
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

//  print all k states
//    for (r = 0; r < kindex; r++)
//    {
//        for(int j = 0; j < col; j++)
//        { cout << arr[r][j];}
//        cout << endl;
//    }
//

    return kindex;


}





