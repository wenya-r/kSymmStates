#include "basicAlg.hpp"
#include "gtest/gtest.h"
#include "rotation.hpp"

#include <vector>

    
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}



TEST(example, sum)
{
    int k = 1;
    ASSERT_EQ(k, 1);
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



TEST(overlapK, kComp)
{
    
    vector<double> state{ -0.02733087971, 0.15430335, -0.1269724703, -0.4355791702,
                          0.2812758202, -0.15430335, -0.1269724703, 0.15430335,
                          0.2812758202, -0.3086066999, 0.2812758202, 0.15430335,
                          -0.1269724703, -0.15430335, 0.2812758202,
                         -0.4355791702, -0.1269724703, 0.15430335, -0.02733087971};

    vector<string> states;
    int k = 0;
    double overlap, numStates_m;
    numStates_m = outputmStates(4, 0, states);
    overlap = overlapKmode(state, k,  19, states) ;
    ASSERT_TRUE(abs(overlap- 0.951634)<0.000001);
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
