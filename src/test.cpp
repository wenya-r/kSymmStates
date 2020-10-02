#include "basicAlg.hpp"
#include "gtest/gtest.h"
#include "rotation.hpp"

#include <vector>
#include <string>
#include <complex>
    
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


TEST(compl3x, mul)
{
    int index;
    complex<double> a(1, 2), b(3, 4), result;
    result = a*b;
    ASSERT_DOUBLE_EQ(result.real(), -5);
}

TEST(compl3x, mulimag)
{
    int index;
    complex<double> a(1, 2), b(3, 4), result;
    result = a*b;
    ASSERT_DOUBLE_EQ(result.imag(), 10);
}

TEST(normalize, sum)
{
    vector<double> state{3.0/5, 4.0/5};

    ASSERT_TRUE(abs(normalize(state, 2)- 1)<0.000000000001);
}


vector<double> site4State()
{
return { -0.02733087971, 0.15430335, -0.1269724703, -0.4355791702,
                          0.2812758202, -0.15430335, -0.1269724703, 0.15430335,
                          0.2812758202, -0.3086066999, 0.2812758202, 0.15430335,
                          -0.1269724703, -0.15430335, 0.2812758202,
                         -0.4355791702, -0.1269724703, 0.15430335, -0.02733087971};
}


TEST(normalizesite4, sum)
{
    vector<double> state;
    state = site4State();
    ASSERT_DOUBLE_EQ(normalize(state,19), 1.000000000059629);
}


TEST(bitfindIndex, length5)
{
    int index;
    vector<int> state{ 64, 60, 40, 13, 4 }; 
    vector<string> vect{"2101", "2020", "1111", "0111", "0011"};
    index = bitToIndex("1111" , vect, state);
    ASSERT_EQ(index, 2);
}


TEST(bitfindIndex, length5_1)
{
    int index;
    vector<int> state{ 64, 60, 40, 13, 4 }; 
    vector<string> vect{"2101", "2020", "1111", "0111", "0011"};
    index = bitToIndex("2101" , vect, state);
    ASSERT_EQ(index, 0);
}


TEST(bitfindIndex, length5_2)
{
    int index;
    vector<int> state{ 64, 60, 40, 13, 4 }; 
    vector<string> vect{"2101", "2020", "1111", "0111", "0011"};
    index = bitToIndex("2020" , vect, state);
    ASSERT_EQ(index, 1);
}


TEST(bitfindIndex, length5_3)
{
    int index;
    vector<int> state{ 64, 60, 40, 13, 4 }; 
    vector<string> vect{"2101", "2020", "1111", "0111", "0011"};
    index = bitToIndex("0111" , vect, state);
    ASSERT_EQ(index, 3);
}


TEST(bitfindIndex, length5_4)
{
    int index;
    vector<int> state{ 64, 60, 40, 13, 4 }; 
    vector<string> vect{"2101", "2020", "1111", "0111", "0011"};
    index = bitToIndex("0011" , vect, state);
    ASSERT_EQ(index, 4);
}



TEST(overlapK, kComp)
{
    
    vector<double> state;
    vector<string> states;
    vector<int> indexTable;
    int k = 0;
    double overlap, numStates_m;
    state = site4State();
    numStates_m = outputmStates(4, 0, states, indexTable);
    
//    cout << "finish numStates_m"<< numStates_m << endl;
    overlap = overlapKmode(state, k, 19, states, indexTable) ;
    ASSERT_DOUBLE_EQ(overlap, 0.95163397553033369);
}

TEST(overlapK, k2Comp)
{
    
    vector<double> state;
    vector<string> states;
    vector<int> indexTable;
    int k = 2;
    double overlap, numStates_m;
    state = site4State();
    numStates_m = outputmStates(4, 0, states, indexTable);
    
//    cout << "finish numStates_m"<< numStates_m << endl;
    overlap = overlapKmode(state, k, 19, states, indexTable) ;
    ASSERT_DOUBLE_EQ(overlap, 0.048366024588924757);
}
TEST(overlapK, imagNorm)
{
    complex<double> *state;
    complex<double> item;
    double norm;
    state = new complex<double>[19]; 
   vector<double> ground;
   ground = site4State();

    for (int i = 0; i< 19; i++)
    { 
        item = {0, ground[i]};
        state[i] = item;
//        cout<< state[i] << endl;
    }
    norm = normalize(state, 19) ;
    ASSERT_DOUBLE_EQ(norm, 1.000000000059629);
}
//
//
//
TEST(overlapK, imagDotProduct)
{
    complex<double> *state;
    complex<double> item;
    cout.precision(17);      double norm;
    state = new complex<double>[19]; 
    vector<double> ground;
    ground = site4State();


    for (int i = 0; i< 19; i++)
    {     state[i] =  {0, ground[i]};    }
    norm = sqrt(dotProduct(state, state, 19)) ;
//    cout << "dotProd : " << dotProduct(state, state, 19)<< endl;
//    cout << " norm :   " <<sqrt(norm) << endl;
//    cout << " normalize" << normalize(state, 19) << endl;
    ASSERT_DOUBLE_EQ(sqrt(norm), 1.000000000059629);
}
//
//
//
//
//
//
//

TEST(overlapK, imagkComp)
{
    complex<double> *state;
    complex<double> item;
    state = new complex<double>[19]; 
    vector<double> ground;
    ground = site4State();                        
    for (int i = 0; i< 19; i++)
    { 
        item = {0, ground[i]};
        state[i] = item;
//        cout<< state[i] << endl;
    }
    vector<string> states;
    vector<int> indexTable;
    int k = 0;
    double overlap, numStates_m;
    numStates_m = outputmStates(4, 0, states, indexTable);
//    cout << "finish numStates_m" << endl;
    overlap = overlapKmode(state, k, 19, states, indexTable) ;
    ASSERT_DOUBLE_EQ(overlap, 0.95163397553033369);
}
//
//

TEST(overlapK, imagk2Comp)
{
    complex<double> *state;
    complex<double> item;
    state = new complex<double>[19]; 
    vector<double> ground;
    ground = site4State();                        
    for (int i = 0; i< 19; i++)
    { 
        item = {0, ground[i]};
        state[i] = item;
//        cout<< state[i] << endl;
    }
    vector<string> states;
    vector<int> indexTable;
    int k = 2;
    double overlap, numStates_m;
    numStates_m = outputmStates(4, 0, states, indexTable);
    cout << "finish numStates_m" << numStates_m << endl;
    overlap = overlapKmode(state, k, 19, states, indexTable) ;
    ASSERT_DOUBLE_EQ(overlap, 0.048366024588924757);
}

//    vector<double> groundState{
//        0,0,0,1/sqrt(6.),0,
//        0,0,1/sqrt(12.),0,1/sqrt(3.0),
//        0,0,0,0,0,
//        1/sqrt(6.), 0,  1/sqrt(12.), 0};


//    vector<double> groundState{
//        0,0,0,0,0,
//        0,0,1/sqrt(6.),0,1/sqrt(3.0),
//        0,0,0,0,0,
//        1/sqrt(3.), 0,  1/sqrt(6.), 0};
