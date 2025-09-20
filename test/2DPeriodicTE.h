#pragma once

#include <vector>
#include <complex>
#include <cmath>


std::vector<double> altComp(std::vector<std::complex<double>>& data);

/*
write a function that a given 1D array of complex values and
total time integration (integer value equal to size of array)
parse into new array twice the length splitting complex
into alternating values of real and imaginary in preparation
for FFT algorithm

then encapsulate FFT algorithm as described and return
array of complex vals for fourier transformed values


*/
void four1(std::vector<double>& data, const int n, const int isign);


/*
write function to find peaks in each FFT of 
*/



/*
write a function that given reference to 2D array of doubles,
writes to csv file for matlab in preparation for pcolor
*/
void writeToFile(std::vector<std::vector<double>>& data);