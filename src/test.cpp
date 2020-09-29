#include "basicAlg.hpp"
#include "gtest/gtest.h"


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

    ASSERT_TRUE(abs(normalize(state, 2)- 1)<0.001);
}

