#include "basicAlg.hpp"
#include "gtest/gtest.h"
#include "rotation.hpp"

#include <vector>
#include <string>
    
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




TEST(normalize, sum)
{
    vector<double> state{3.0/5, 4.0/5};

    ASSERT_TRUE(abs(normalize(state, 2)- 1)<0.000000000001);
}

TEST(normalizesite4, sum)
{
    vector<double> state{ -0.02733087971, 0.15430335, -0.1269724703, -0.4355791702,
                          0.2812758202, -0.15430335, -0.1269724703, 0.15430335,
                          0.2812758202, -0.3086066999, 0.2812758202, 0.15430335,
                          -0.1269724703, -0.15430335, 0.2812758202,
                         -0.4355791702, -0.1269724703, 0.15430335, -0.02733087971};
    ASSERT_TRUE(abs(normalize(state, 19)- 1)<0.000000001);
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
    
    vector<double> state{ -0.02733087971, 0.15430335, -0.1269724703, -0.4355791702,
                          0.2812758202, -0.15430335, -0.1269724703, 0.15430335,
                          0.2812758202, -0.3086066999, 0.2812758202, 0.15430335,
                          -0.1269724703, -0.15430335, 0.2812758202,
                         -0.4355791702, -0.1269724703, 0.15430335, -0.02733087971};

    vector<string> states;
    vector<int> indexTable;
    int k = 0;
    double overlap, numStates_m;
    numStates_m = outputmStates(4, 0, states, indexTable);
//    cout << "finish numStates_m" << endl;
    overlap = overlapKmode(state, k, 19, states, indexTable) ;
    ASSERT_TRUE(abs(overlap- 0.951634)<0.000001);
}

//
//
//
//TEST(bitFindIndex, length6)
//{
//    vector<double> state{ -0.02733087971, 0.15430335, -0.1269724703, -0.4355791702,
//                          0.2812758202, -0.15430335, -0.1269724703, 0.15430335,
//                          0.2812758202, -0.3086066999, 0.2812758202, 0.15430335,
//                          -0.1269724703, -0.15430335, 0.2812758202,
//                         -0.4355791702, -0.1269724703, 0.15430335, -0.02733087971};
//    ASSERT_TRUE(abs(normalize(state, 19)- 1)<0.000000001);
//}
//
//
//
//TEST(bitFindIndex, length7)
//{
//    vector<double> state{ -0.02733087971, 0.15430335, -0.1269724703, -0.4355791702,
//                          0.2812758202, -0.15430335, -0.1269724703, 0.15430335,
//                          0.2812758202, -0.3086066999, 0.2812758202, 0.15430335,
//                          -0.1269724703, -0.15430335, 0.2812758202,
//                         -0.4355791702, -0.1269724703, 0.15430335, -0.02733087971};
//    ASSERT_TRUE(abs(normalize(state, 19)- 1)<0.000000001);
//}


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
