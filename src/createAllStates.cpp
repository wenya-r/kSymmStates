#include <iostream>
#include <cmath>
#include <string>

using namespace std;

void outputAllStates(int L);
void outputmStates(int L, int m);


int main()
{
    int L = 4;
   outputmStates(L, 0);

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

void outputmStates(int L, int m){
// m is Sz total. sum(string)= L+m 

    int numOfStates = pow(3,L);
    int stateNum, bit, sumBit;
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
            
            cout << state << endl;
        }
    }


}



