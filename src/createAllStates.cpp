#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>

using namespace std;

void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect);
int bitToNum(string bit);
string translation(string bit);

double dotProduct(complex<double> * vec1, complex<double> * vec2, int col);
double dotProduct(complex<double> * vec1, vector<double> * vec2, int col);
int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect);

int main()
{
    int L = 4, kindex;
    int num, numStates_m; 
    double overlap = 0, dot = 0;
    vector<string> states;
    numStates_m = outputmStates(L, 0, states);


    complex<double> **A = new complex<double>*[numStates_m];
    for(int i = 0; i < numStates_m; i++)
    {
        A[i] = new complex<double>[numStates_m];
    }
    for (int i = 0; i < numStates_m; i++)
    {
        for (int j = 0; j < numStates_m; j++) {A[i][j] = 0;}
    }

    vector<double> groundState{-0.07453559924999242, 0.22360679774997777, -0.14907119849998562,
-0.44721359549995865, 0.2236067977499792, -0.07453559924999299, -0.14907119849998546,
0.22360679774997771, 0.2236067977499792, -0.2981423969999716, 0.22360679774997913,
0.2236067977499794, -0.1490711984999861, -0.07453559924999298, 0.22360679774997905,
-0.44721359549995787, -0.1490711984999861, 0.2236067977499792, -0.07453559924999312};

    kindex = coefficientskMomentum(A, numStates_m, numStates_m, L, 0, states);
    cout << "kinde: " << kindex << endl;
    for (int i = 0; i < kindex; i++)
    {
        dot = dotProduct(A[i], groundState, numStates_m);
        overlap = overlap + dot*dot;
    }
    cout << "overlap : " << overlap << endl;

    return 0;
} //  end Main




void outputAllStates(int L){
    int numOfStates = pow(3,L);
    int stateNum;
    string state;
    cout << "num of states: " << numOfStates << endl;
    for(int i = numOfStates-1;i>-1;i--)
    {   
        cout << i;
        stateNum = i;
        state = to_string(stateNum%3);
        for (int j = 0; j < 3; j++)
        {
            stateNum = stateNum/3;
            state = to_string(stateNum%3) + state;
        }
        cout << state << endl;
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
        for (int j = 0; j < 3; j++)
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
    cout << "do you see bit?" << bit<< endl;
    L = bit.length();
    num = stoi(bit, nullptr, 3);
//    for (int i = 1; i < L  ; i++)
//    {
//        num = num*3;
//        num = num + stoi(bit[i]) ;
//    }
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
//    cout << " bit is " << bit << endl;
//    cout << transBit << endl;
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


double dotProduct(complex<double> * vec1, vector<double> * vec2, int col)
{
    complex<double> product = 0;
    for (int i = 0; i < col; i++)
    {        product = product + vec1[i].real()*(vec2[i]) + i*vec1[i].imag()*(vec2[i]); }
    return norm(product);
}





int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect)
{
    string bit;
    int i, r, kindex = 0;
    double km, norm, product ;
    bool iskState;
    const double pi = acos(-1);
    complex<double> * cofArr;
  
    //declare complex array and set to zero
    cofArr = new complex<double>[col];
    for (int i = 0; i < col; i++) {cofArr[i]=0;}
    cout << " the number of cofAff = " << sizeof(cofArr) << endl;
//
    km = k*2*pi/vect[0].length();
    
    for(i = 0; i < row; i++)
    {
        cout << i << endl;
        bit = vect[i];
        for (int i = 0; i < col; i++) {cofArr[i]=0;}
        //construct the first state
        // should be a function find out the coeff for |a> and k 
        
        
        cofArr[bitToIndex(bit, vect)] = 1;        
        for(r = 1; r < L; r++)
        {
            bit = translation(bit);
           
            cofArr[bitToIndex(bit, vect)] =  cofArr[bitToIndex(bit, vect)]+exp(i*km*r);        
        }

//        for(r = 0; r < col ; r++)
//        {
//            cout << r << cofArr[r] ;
//        }
//        cout << endl;         
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
    for (r = 0; r < kindex; r++)
    {
        for(int j = 0; j < col; j++)
        { cout << arr[r][j];}
        cout << endl;
    }


    return kindex;




//    int i, j, k;
//    for (i = 0 ; i < row; i++){
//        for (j = 0; j < col; j++){
//            arr[i][j] = 8+j;
//        }
//    }  

}





