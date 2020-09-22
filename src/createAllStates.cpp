#include <iostream>
#include <cmath>
#include <string>

using namespace std;

void outputAllStates(int L);
int outputmStates(int L, int m);
int bitToNum(string bit);
string translation(string bit);

void coefficientskMomentum(double coef[][10]);

int main()
{
    int L = 4;
    int num;
    double coeffs[19][10];
//    double & coeffsRef = coeffs;
    cout << "row" << sizeof(coeffs)/sizeof(coeffs[0]) << endl;
    cout << "column" << sizeof(coeffs[0])/sizeof(coeffs[0][0]) << endl;
    
    coefficientskMomentum(coeffs);
    cout << "coeffs : " << coeffs[5][6] << endl;
//    num = outputmStates(L, 0);
//   num = bitToNum("0022");
    cout << " num of states: "<< num << endl;
    return 0;
}

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

int outputmStates(int L, int m){
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
            cout << state << endl;
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
    return "0";

}


void coefficientskMomentum(double coef[][10])
{
    
//    cout << "row" << sizeof(coeffs)/sizeof(coeffs[0]) << endl;
//    cout << "column" << sizeof(coeffs[0])/sizeof(coeffs[0][0]) << endl;
    int row, col;
    row = sizeof(coef)/sizeof(coef[0]);
    col = sizeof(coef[0])/sizeof(coef[0][0]) ;
    for (int i = 0 ; i < row; i++){
        for (int j = 0; j < col; j++){
            coef[i][j] = 8;
        }
    }  

}






