#include "basicAlg.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <complex>
#include <fstream>


using namespace std;



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




