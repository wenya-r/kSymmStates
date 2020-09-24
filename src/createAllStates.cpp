#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include "H5Cpp.h"

#define FILE1 "highSpin_6sites_m0_groundState.hdf5"
#define DATASET1 "E"
#define DATASET0 "L"

#define size_vec 3241135

using namespace std;

void outputAllStates(int L);
int outputmStates(int L, int m, vector<string> &vect);
int bitToNum(string bit);
string translation(string bit);

double dotProduct(complex<double> * vec1, complex<double> * vec2, int col);
double dotProduct(complex<double> * vec1, vector<double> vec2, int col);
int coefficientskMomentum(complex<double> **arr, int row, int col, int L, int k, vector<string> &vect);


const H5std_string FILE_NAME("highSpin_6sites_m0_groundState.hdf5");
const H5std_string DATASET1_NAME("E");
const H5std_string DATASET2_NAME("L");
const int dim_num = 1000;


int main()
{

    int L = 4, kindex;
    int num, numStates_m; 
    double overlap = 0, dot = 0;
    vector<string> states;
    numStates_m = outputmStates(L, 0, states);

    H5::H5File file( FILE_NAME, H5F_ACC_RDONLY );
    H5::DataSet dataset1 = file.openDataSet( DATASET1_NAME ); 
    H5::DataSet dataset2 = file.openDataSet( DATASET2_NAME ); 

    H5T_class_t type_class1 = dataset1.getTypeClass();
    H5T_class_t type_class2 = dataset2.getTypeClass();

    H5::IntType intype = dataset1.getIntType();




    complex<double> **A = new complex<double>*[numStates_m];
    for(int i = 0; i < numStates_m; i++)
    {
        A[i] = new complex<double>[numStates_m];
    }
    for (int i = 0; i < numStates_m; i++)
    {
        for (int j = 0; j < numStates_m; j++) {A[i][j] = 0;}
    }

    vector<double> groundState{
-0.027330879711206537,
0.1543033499620927,
-0.12697247025088626,
-0.43557917017507153,
0.2812758202129786,
-0.15430334996209233,
-0.1269724702508861,
0.15430334996209158,
0.2812758202129785,
-0.3086066999241836,
0.28127582021297703,
0.1543033499620921,
-0.1269724702508852,
-0.15430334996209222,
0.28127582021297703,
-0.43557917017506625,
-0.1269724702508852,
0.1543033499620907,
-0.027330879711206155

     };


//    vector<double> groundState{-0.07453559924999242, 0.22360679774997777, -0.14907119849998562,
//-0.44721359549995865, 0.2236067977499792, -0.07453559924999299, -0.14907119849998546,
//0.22360679774997771, 0.2236067977499792, -0.2981423969999716, 0.22360679774997913,
//0.2236067977499794, -0.1490711984999861, -0.07453559924999298, 0.22360679774997905,
//-0.44721359549995787, -0.1490711984999861, 0.2236067977499792, -0.07453559924999312};
   
//    vector<double> groundState{
//0.0030355420235822026,
//-0.016114817125276463,
//0.02098523773539049,
//-0.007905962633696056,
//0.03894193584436984,
//-0.04381235645448381,
//0.02098523773539064,
//0.03894193584436937,
//-0.016114817125276605,
//0.0030355420235820586,
//0.02098523773539046,
//-0.04500601749436345,
//0.020985237735390436,
//-0.043812356454483864,
//0.10694145627338568,
//-0.06799952042901564,
//-0.1539557389670103,
//0.07103506245259766,
//-0.01611481712527653,
//-0.08915814477713649,
//0.09402856538725017,
//0.09402856538725012,
//-0.06799952042901528,
//0.02098523773539004,
//0.020985237735390003,
//-0.007905962633695938,
//0.03894193584436971,
//-0.15395573896701034,
//0.09402856538725056,
//0.3529174954283832,
//-0.15395573896701037,
//0.0389419358443697,
//0.09402856538725012,
//-0.08915814477713625,
//-0.15395573896700998,
//0.10694145627338451,
//-0.043812356454483364,
//-0.04500601749436228,
//0.02098523773538999,
//0.03894193584436947,
//-0.043812356454483226,
//0.038941935844369197,
//0.02098523773538992,
//-0.016114817125276154,
//0.003035542023582111,
//-0.007905962633696039,
//0.02098523773539045,
//-0.016114817125276456,
//0.02098523773539048,
//-0.06799952042901565,
//0.07103506245259773,
//0.09402856538725055,
//-0.06799952042901564,
//0.020985237735390475,
//0.09402856538725018,
//-0.15395573896701018,
//-0.08915814477713625,
//0.10694145627338489,
//-0.04500601749436266,
//-0.043812356454483385,
//0.020985237735390048,
//-0.016114817125276536,
//0.07103506245259766,
//-0.06799952042901565,
//-0.15395573896701034,
//0.10694145627338572,
//-0.043812356454483864,
//-0.06799952042901528,
//0.10694145627338489,
//0.10694145627338453,
//-0.14891893414133686,
//0.10694145627338486,
//0.10694145627338451,
//-0.06799952042901526,
//-0.043812356454483226,
//0.10694145627338422,
//-0.15395573896700887,
//-0.06799952042901497,
//0.0710350624525968,
//-0.016114817125276154,
//0.020985237735389996,
//-0.04381235645448338,
//-0.04500601749436229,
//0.10694145627338453,
//-0.08915814477713624,
//-0.15395573896700998,
//0.09402856538725012,
//0.02098523773538992,
//-0.06799952042901497,
//0.09402856538724944,
//0.07103506245259708,
//-0.06799952042901496,
//0.02098523773538992,
//-0.01611481712527626,
//0.02098523773538998,
//-0.007905962633695888,
//0.003035542023582107,
//-0.016114817125276536,
//0.02098523773539048,
//0.038941935844369696,
//-0.04381235645448385,
//0.038941935844369904,
//0.020985237735390048,
//-0.04500601749436266,
//-0.04381235645448338,
//0.10694145627338487,
//-0.15395573896701015,
//-0.08915814477713624,
//0.09402856538725016,
//0.0389419358443692,
//-0.1539557389670089,
//0.35291749542837997,
//0.09402856538724944,
//-0.15395573896700887,
//0.038941935844369197,
//-0.007905962633695938,
//0.020985237735390048,
//0.020985237735389996,
//-0.06799952042901528,
//0.09402856538725017,
//0.0940285653872501,
//-0.08915814477713649,
//-0.016114817125276154,
//0.0710350624525968,
//-0.15395573896700887,
//-0.06799952042901497,
//0.10694145627338422,
//-0.043812356454483226,
//0.02098523773538998,
//-0.045006017494362276,
//0.020985237735389975,
//0.0030355420235821085,
//-0.016114817125276154,
//0.038941935844369197,
//0.02098523773538992,
//-0.043812356454483226,
//0.03894193584436948,
//-0.007905962633695886,
//0.02098523773538998,
//-0.01611481712527627,
//0.0030355420235821454    };
    for (int i = 0; i < numStates_m; i++)
    {    cout << groundState[i] << endl;}
    cout << endl;
    kindex = coefficientskMomentum(A, numStates_m, numStates_m, L, 0, states);
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
        cout << i << endl;
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





