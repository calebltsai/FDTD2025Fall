#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>

#include "2DPeriodicTE.h"
#include "2DPeriodicTE.cpp"

// start by defining constants

double pi = 3.14159265358979323846;

double c = 299792458;
double mu = pi*4.0e-7;
double epsilon = 1.0/mu/c/c;

unsigned resolution = 100;
double courantfactor = 0.5;
double freq = 1.0e9;
double omega = 2.0*pi*freq;
double lambda = 2.0*pi*c/omega;

double dx = lambda/resolution;
double dt = courantfactor*dx/c;
double dy = dx;

// define grid boundaries and simulation parameters
const int ib = 64;
const int jb = 64;
const int totaltime = std::pow(2,10+4);

std::complex<double> Ez[ib+1][jb+1];
std::complex<double> Hy[ib+1][jb+1];
std::complex<double> Hx[ib+1][jb+1];

double eps[ib+1][ib+1];
double mus[ib+1][ib+1];

// index of refraction for n1 outside, n2 inside
unsigned n1 = 1;
unsigned n2 = 2;


// setting brillouin zone and traversal
const int range = 500;
const int numK = 10;

double bandstructure[numK][range];

int main() {   
	// setting dielectrics
	for (int x=0;x<=ib+1;x++){
		for (int y=0;y<=jb+1;y++){
			if ( (x-ib/2)*(x-ib/2)+(y-ib/2)*(y-ib/2) >= (ib/4)*(ib/4) ) eps[x][y] = epsilon*n2*n2;
			else eps[x][y] = epsilon*n1*n1;
		}
	}
    
 
    for (int k=0; k < 3*(numK+1); k++){
		for (int t=0; t < totaltime ;t++){
			for (int x=0; x < ib; x++) {
				for (int y=0; y < jb; y++) {
					
				}
			}
		}
		std::cout << k << std::endl;
	}

/*
% constants

% handle the array of k vectors in reduced brillouin zone used for bloch condition

kxs = [linspace(0, pi/(dx*ib), numK+1), ones(1, numK+1)*pi/(dx*ib), linspace(pi/(dy*ib), 0, numK+1)]' ;
kys = [zeros(1,numK+1), linspace(0, pi/(dy*ib), numK+1), linspace(pi/(dy*ib), 0, numK+1)]';

% core time evolution simulation
for k=1:(numK+1)*3
    k
    e_obs = zeros(totaltime, 1);
    % time step evolution
    for time=1:totaltime
        % update e with previous e data and 1/2 frame data of h
        for x=2:ib+1
            for y=2:jb+1
                HxDy = ( Hx(x,y)-Hx(x,y-1) )/dy;
                HyDx = ( Hy(x,y)-Hy(x-1,y) )/dx;
    
                Ez(x,y) = Ez(x,y) + dt/epss(x,y)*( HyDx - HxDy );
            end
        end
        Ez(1,:)=Ez(ib+1,:)*exp(1i*kxs(k, 1)*ib*dx);
        Ez(:,1)=Ez(:,jb+1)*exp(1j*kys(k, 1)*jb*dy);

        Ez(ib/2, jb/2) = exp(-(time-10)^2/((2)^2/10))*sin(omega*dt*(time-10)) + Ez(ib/2,jb/2);
    
        % scattered boundary correction
        % update h with previous h data and 1/2 frame data of e

        for x=1:ib
            for y=1:jb
                % total/scattered fields
                EzDy = ( Ez(x,y+1)-Ez(x,y) )/dy;
                EzDx = ( Ez(x+1,y)-Ez(x,y) )/dx;
    
                Hy(x,y) = Hy(x,y) + dt/mus(x,y)*( EzDx );
                Hx(x,y) = Hx(x,y) - dt/mus(x,y)*( EzDy );
            end
        end
        Hx(:,jb+1) = Hx(:,1)*exp(-1j*kys(k,1)*jb*dy);
        Hy(ib+1,:) = Hy(1,:)*exp(-1i*kxs(k,1)*ib*dx);

        e_obs(time, 1) = Ez(ib/4, jb/4);
    end

    transform = fft(e_obs, totaltime)/totaltime;
    f=1/2/dt*linspace(0,1,totaltime/2+1)';

    bandstructure(k,1:range) = abs(transform(1:range,1));
end

*/
	
	return 0;
}