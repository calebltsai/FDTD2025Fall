#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>

#include "2DPeriodicTM.h"
#include "2DPeriodicTM.cpp"

// CHANGES TO MAKE:
/* 
optional:
1. reduce use of push_back
2. use arrays instead of vectors
3. allocate all the memory needed in the beginning so that
there is no need for implicit resizing operations, esp in vectors
4. use gdb and valgrind to catch and throw errors

practice using heap memory for temporary variables?
does anything make more sense for dynamic memory?
*/
// start by defining constants
const double pi = 3.14159265358979323846;

double c = 299792458;
double mu = pi*4.0e-7;
double epsilon = 1.0/mu/c/c;

double resolution = 100;
double courantfactor = 0.5;
double freq = 1.0e9;
double omega = 2.0*pi*freq;
double lambda = 2.0*pi*c/omega;

double dx = lambda/resolution;
double dt = courantfactor*dx/c;
double dy = dx;

// define grid boundaries and simulation parameters
const int ib = 32;
const int jb = 32;

// initialize to complex value zero
std::complex<double> Hz[ib+1][jb+1] = {std::complex<double>(0,0)};
std::complex<double> Ey[ib+1][jb+1] = {std::complex<double>(0,0)};
std::complex<double> Ex[ib+1][jb+1] = {std::complex<double>(0,0)};

double eps[ib+1][jb+1];
double mus[ib+1][jb+1];

// index of refraction for n1 outside, n2 inside
int n1 = 1;
int n2 = 2;

// setting brillouin zone and traversal
const int numK = 50;
const int tmp = std::pow(2,10+4);

std::vector<std::vector<double>> bandstructure(numK*3, std::vector<double>(tmp));

int main() {   
	// setting dielectrics and dimagnetics
	for (int x=0;x<=ib;x++){
		for (int y=0;y<=jb;y++){
			mus[x][y] = mu;
			// normal/vacuum region inside (cavity)
			if ( (x-ib/2)*(x-ib/2)+(y-ib/2)*(y-ib/2) >= (ib/3)*(ib/3) ) eps[x][y] = epsilon*n2*n2;
			else eps[x][y] = epsilon*n1*n1;
		}
	}

	std::ofstream outfile("test.csv");

	// begin analysis by iteration over reduced brillouin zone
 double kxs[numK*3];
	double kys[numK*3];

	for (int idx=0; idx<numK; idx++) {
		kxs[idx] = pi/(ib*dx)*idx/numK;
		kys[idx] = 0;
		kxs[idx+numK+1] = pi/(ib*dx);
		kys[idx+numK+1] = pi/(jb*dy)*idx/numK;
		kxs[idx+2*(numK+1)] = pi/(ib*dx)*(numK-idx)/numK;
		kys[idx+2*(numK+1)] = pi/(jb*dy)*(numK-idx)/numK;
	}


	for (int k=0; k < 3*numK; k++) {
		std::cout << "k=" << k << std::endl;

		const int totaltime = std::pow(2,10+4);
		std::vector<std::complex<double>> hobs(totaltime);

		std::complex<double> Hz[ib+1][jb+1] = {std::complex<double>(0,0)};
		std::complex<double> Ey[ib+1][jb+1] = {std::complex<double>(0,0)};
		std::complex<double> Ex[ib+1][jb+1] = {std::complex<double>(0,0)};

		for (int t=0; t < totaltime; t++) {
			
			// Ex, Ey
			for (int x=0; x <= ib-1; x++) {
				for (int y=0; y <= jb-1; y++) {
					std::complex<double> HzDy = ( Hz[x][y+1]-Hz[x][y] )/dy;
     std::complex<double> HzDx = ( Hz[x+1][y]-Hz[x][y] )/dx;

					Ex[x][y] = Ex[x][y] + dt/eps[x][y]*( HzDy );
					Ey[x][y] = Ey[x][y] - dt/eps[x][y]*( HzDx );
				}
			}
			// loop over boundary edge for bloch periodic boundary condition
			for (int idx = 0;idx<ib+1;idx++) {
				Ey[ib][idx] = Ey[0][idx]*std::exp(std::complex<double>(0,-1)*kxs[k]*std::complex<double>(ib,0)*dx);
				Ex[idx][jb] = Ex[idx][0]*std::exp(std::complex<double>(0,-1)*kys[k]*std::complex<double>(jb,0)*dy);
			}

			// source pulse not from cavity
			// Hz
			for (int x=1; x <= ib; x++) {
				for (int y=1; y <= jb; y++) {
					std::complex<double> ExDy = ( Ex[x][y]-Ex[x][y-1] )/dy;
     std::complex<double> EyDx = ( Ey[x][y]-Ey[x-1][y] )/dx;

					Hz[x][y] = Hz[x][y] - dt/mus[x][y]*( EyDx - ExDy );
				}
			}
			
			// bloch boundary condition on Hz
			for (int idx = 0;idx<ib+1;idx++) {
				Hz[0][idx]=Hz[ib][idx]*std::exp(std::complex<double>(0,1)*kxs[k]*std::complex<double>(ib,0)*dx);
				Hz[idx][0]=Hz[idx][jb]*std::exp(std::complex<double>(0,1)*kys[k]*std::complex<double>(jb,0)*dy);
			}
			
			Hz[ib/2+5][jb/2+5] = std::exp(-std::pow((t-100),2)/std::pow(20,2)/10)*std::exp(std::complex<double>(0,omega*dt*(t-100))) + Hz[ib/2+5][jb/2+5];

			hobs[t] = Hz[ib/2-5][jb/2-5];
			if (outfile.is_open()){
				outfile << hobs[t] << ", t=" << t << "\n";
			}		
		}

		std::vector<double> transform = altComp(hobs);
		four1(transform, totaltime, 1);

		// function: take abs val of tuple complex coords in the transform, then pass back to bandstructure
		std::vector<double> temp;
		for (int idx = 0;idx<totaltime*2;idx+=2) {
			// taking the norm of each pair of values in the array
			temp.push_back(std::sqrt(transform.at(idx)*transform.at(idx)+transform.at(idx+1)*transform.at(idx+1)));
		}
  
		bandstructure[k] = temp;
	}

	writeToFile(bandstructure);

	return 0;
}