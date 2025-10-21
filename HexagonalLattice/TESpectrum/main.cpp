#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>

#include "2DPeriodicTE.h"
#include "2DPeriodicTE.cpp"

/* 
practice using heap memory for temporary variables? (does anything make more sense as dynamic memory?)

splitting into functions (and classes?) for organizing:
complex operations?
matrix operations?
time step
setting boundaries/pml/tfsf

valgrind (using and managing memory through pointers, etc)
pointer of array<complex> with dimension of rows->pointing to column element-> which is a pointer to complex value
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
// jb = 54*sqrt(3)
const int ib = 54;
const int jb = 94;
const int a0 = ib/3;
const double arRatio[3] = {3.0, 2.83, 3.17};

// center points: middle, tl, bl, tr, br
const int ic[] = {ib/2, 1, 1, ib, ib};
const int jc[] = {jb/2, 1, jb, 1, jb};

// initialize to complex value zero
std::complex<double> Ez[ib+1][jb+1] = {std::complex<double>(0,0)};
std::complex<double> Hy[ib+1][jb+1] = {std::complex<double>(0,0)};
std::complex<double> Hx[ib+1][jb+1] = {std::complex<double>(0,0)};

std::vector<std::vector<double>> eps(ib+1, std::vector<double>(jb+1,epsilon));
std::vector<std::vector<double>> mus(ib+1, std::vector<double>(jb+1,mu));

// index of refraction for n1 outside, n2 inside
int n1 = 1;
int n2 = 2;

// setting brillouin zone and traversal
const int numK = 40;
const int tmp = std::pow(2,10+4);

std::vector<std::vector<double>> bandstructure(numK*3, std::vector<double>(tmp));

int main() {
	// geometry of the cylindrical dielectrics with respect to a center (0,0)
	double increment = 2.0*pi/6.0;
	int a[6][2] = {
		{static_cast<int>(std::round(a0*std::cos(increment*0.0))), static_cast<int>(std::round(a0*std::sin(increment*0.0)))},
		{static_cast<int>(std::round(a0*std::cos(increment*1.0))), static_cast<int>(std::round(a0*std::sin(increment*1.0)))},
		{static_cast<int>(std::round(a0*std::cos(increment*2.0))), static_cast<int>(std::round(a0*std::sin(increment*2.0)))},
		{static_cast<int>(std::round(a0*std::cos(increment*3.0))), static_cast<int>(std::round(a0*std::sin(increment*3.0)))},
		{static_cast<int>(std::round(a0*std::cos(increment*4.0))), static_cast<int>(std::round(a0*std::sin(increment*4.0)))},
		{static_cast<int>(std::round(a0*std::cos(increment*5.0))), static_cast<int>(std::round(a0*std::sin(increment*5.0)))}
	};

	// setting dielectrics and dimagnetics
	for (int x=0;x<=ib;x++){
		for (int y=0;y<=jb;y++){
			// mus[x][y] = mu;
			// normal/vacuum region outside (scattering obstacle)
			for (int i = 0; i < 5; i++) {
				for (int j = 0; j < 6; j++) {
					if ( ((ic[i]+a[j][0])-x)*((ic[i]+a[j][0])-x)+((jc[i]+a[j][1])-y)*((jc[i]+a[j][1])-y) < a0*a0/(arRatio[0]*arRatio[0]) ) eps[x][y] = epsilon*n2*n2;
				}
			}
		}
	}
	writeToFile(eps, "output/dielectrics.csv");

	// begin analysis by iteration over reduced brillouin zone
 double kxs[numK*3];
	double kys[numK*3];

	for (int idx=0; idx<numK; idx++) {
		kxs[idx] = 4.0*pi/(3.0*ib*dx)*(idx*1.0/numK*1.0);
		kys[idx] = 0;
		kxs[idx+numK] = 4.0*pi/(3.0*ib*dx)*(1.0-(idx*1.0/numK*1.0)/4.0);
		kys[idx+numK] = pi/(std::sqrt(3)*ib*dy)*idx/(1.0*numK);
		kxs[idx+2*numK] = pi/(ib*dx)*(numK-idx)/(1.0*numK);
		kys[idx+2*numK] = pi/(std::sqrt(3)*ib*dy)*(numK-idx)/(numK*1.0);
	}

	for (int k=0; k < 3*numK; k++) {
		std::cout << "k=" << k << std::endl;

		const int totaltime = std::pow(2,10+4);
		std::vector<std::complex<double>> eobs(totaltime);

		std::complex<double> Hz[ib+1][jb+1] = {std::complex<double>(0,0)};
		std::complex<double> Ey[ib+1][jb+1] = {std::complex<double>(0,0)};
		std::complex<double> Ex[ib+1][jb+1] = {std::complex<double>(0,0)};

		for (int t=0; t < totaltime; t++) {
			
			// Ez
			for (int x=1; x <= ib; x++) {
				for (int y=1; y <= jb; y++) {
					std::complex<double> HxDy = ( Hx[x][y]-Hx[x][y-1] )/dy;
     std::complex<double> HyDx = ( Hy[x][y]-Hy[x-1][y] )/dx;

					Ez[x][y] = Ez[x][y] + dt/eps[x][y]*( HyDx - HxDy );
				}
			}
			// loop over boundary edge for bloch periodic boundary condition
			for (int idx = 0;idx<ib+1;idx++) {
				Ez[0][idx] = Ez[ib][idx]*std::exp(std::complex<double>(0,1)*kxs[k]*std::complex<double>(ib,0)*dx);
				Ez[idx][0] = Ez[idx][jb]*std::exp(std::complex<double>(0,1)*kys[k]*std::complex<double>(jb,0)*dy);
			}

			// source pulse from center
			Ez[ic[0]+3][jc[0]+3] = std::exp(-std::pow((t-100),2)/std::pow(20,2)/10)*std::exp(std::complex<double>(0,omega*dt*(t-100))) + Ez[ic[0]+3][jc[0]+3];
			Ez[ic[1]+3][jc[1]+3] = std::exp(-std::pow((t-100),2)/std::pow(20,2)/10)*std::exp(std::complex<double>(0,omega*dt*(t-100-kxs[k]*ib*dx/2-kys[k]*jb*dy/2))) + Ez[ic[1]+3][jc[1]+3];
			// Hy, Hx
			for (int x=0; x <= ib-1; x++) {
				for (int y=0; y <= jb-1; y++) {
					std::complex<double> EzDy = ( Ez[x][y+1]-Ez[x][y] )/dy;
     std::complex<double> EzDx = ( Ez[x+1][y]-Ez[x][y] )/dx;

					Hy[x][y] = Hy[x][y] + dt/mus[x][y]*( EzDx );
     Hx[x][y] = Hx[x][y] - dt/mus[x][y]*( EzDy );
				}
			}
			
			// bloch boundary condition on Hy, Hx
			for (int idx = 0;idx<ib+1;idx++) {
				Hy[ib][idx]=Hy[0][idx]*std::exp(std::complex<double>(0,-1)*kxs[k]*std::complex<double>(ib,0)*dx);
				Hx[idx][jb]=Hx[idx][0]*std::exp(std::complex<double>(0,-1)*kys[k]*std::complex<double>(jb,0)*dy);
			}

			eobs[t] = Ez[ib/4][jb/4];
		}

		std::vector<double> transform = altComp(eobs);
		four1(transform, totaltime, 1);

		// function: take abs val of tuple complex coords in the transform, then pass back to bandstructure
		std::vector<double> temp;
		for (int idx = 0;idx<totaltime*2;idx+=2) {
			// taking the norm of each pair of values in the array
			temp.push_back(std::sqrt(transform.at(idx)*transform.at(idx)+transform.at(idx+1)*transform.at(idx+1)));
		}
	
		bandstructure[k] = temp;
	}

	writeToFile(bandstructure, "output/2DTESpectrum.csv");
	return 0;
}