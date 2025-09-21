#include "2DPeriodicTE.h"

#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iostream>

const double PI = 3.14159265358979323846;

std::vector<double> altComp(std::vector<std::complex<double>>& data) {
	std::vector<double> ret;
	for (unsigned idx=0;idx< data.size();idx++) {
		ret.push_back(data.at(idx).real());
		ret.push_back(data.at(idx).imag());
	}
	return ret;
}

/*
write a function that a given 1D array of complex values and
total time integration (integer value equal to size of array)
parse into new array twice the length splitting complex
into alternating values of real and imaginary in preparation
for FFT algorithm

then encapsulate FFT algorithm as described and return
array of complex vals for fourier transformed values
*/
void four1(std::vector<double>& data, const int n, const int isign) {
	int nn,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) throw("n must be power of 2 in four1");
	nn = n << 1;
	j = 1;
	for (i=1;i<nn;i+=2) {
		if (j > i) {
			std::swap(data[j-1],data[i-1]);
			std::swap(data[j],data[i]);
		}
		m=n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}


/*
write function to find peaks in each FFT, in the case youare not returning for pcolor
*/



/*
write a function that given reference to 2D array of doubles,
writes to csv file for matlab in preparation for pcolor
*/
void writeToFile(std::vector<std::vector<double>>& data) {
	std::ofstream outFile("output/2DTESpectrum.csv");
	if (outFile.is_open()) {
		for (unsigned i=0; i < data.size();i++) {
			for (unsigned j=0; j < data[0].size();j++) {
				outFile << data[i][j];
				if (j < data[i].size() - 1) {
					outFile << ", ";
				}
			}
			outFile << std::endl;
		}
		outFile.close();
		std::cout << "Sucessful write to file 2DTESpectrum.csv" << std::endl;
	} else std::cout << "Error in opening file 2DTESpectrum.csv" << std::endl;
}