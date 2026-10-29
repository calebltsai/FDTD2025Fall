clear;

c = 299792458;
mu = pi*4.0e-7;
epsilon = 1/mu/c/c;

resolution = 100;
courantfactor = 0.5;
freq = 1e9;
omega = 2*pi*freq;
lambda = 2*pi*c/omega;

dx = lambda/resolution;
dt = courantfactor*dx/c;
totaltime = (2^10)*(2^4);

% portion lengths in lattice geometry
La = 10;
Lb = 40;
ib = La + Lb;
% indices of refraction in lattice geometry
Na = 4;
Nb = 1;

% set k iterations
NumK = 10;
kx = linspace(0, pi/ib/dx, NumK+1)';

% set coefficients
eps = ones(ib, 1)*epsilon;
mus = ones(ib, 1)*mu;

% setting up photonic properties of unit lattice block
eps(1:La, 1) = epsilon*Na^2;
eps(La+1:Lb, 1) = epsilon*Nb^2;

Ca = ones(ib, 1);
Da = ones(ib, 1);
Cb = ones(ib, 1);
Db = ones(ib, 1);
for i=1:ib
    Cb(i, 1) = dt/eps(i, 1);
    Db(i, 1) = dt/mus(i, 1);
end

range = 500;

bandstructure = zeros(NumK, range);
% iteration over ks in the brillouin zone to find highest density of
% states/edge state associated wtih resonant frequency to identify photonic band structure
for k=1:NumK
    % set grids
    Ez = zeros(ib, 1);
    Hy = zeros(ib, 1);
    Eobs = zeros(totaltime,1);
    obs = La+Lb/2;
    for t=1:totaltime
        Ez(1, 1) = Cb(1, 1)/dx*( Hy(1, 1) - Hy(ib, 1)*exp(-1i*kx(k, 1)*ib*dx) )+Ca(1, 1)*Ez(1, 1);
        for i=2:ib
            Ez(i,1) = Cb(i, 1)/dx*( Hy(i,1) - Hy(i-1,1) ) + Ca(i, 1)*Ez(i,1);
        end
        Ez(ib/2, 1) = exp(-(t-10)^2/((2)^2/10))*sin(omega*dt*(t-10)) + Ez(ib/2,1);
        
        %update h with previous h data and 1/2 frame data of e
        for i=1:ib-1
            Hy(i,1) = Db(i, 1)/dx*( Ez(i+1,1) - Ez(i,1) )+Da(i, 1)*Hy(i,1); 
        end
        Hy(ib, 1) = Db(1, 1)/dx*( Ez(1, 1)*exp(1i*kx(k,1)*ib*dx) - Ez(ib, 1) )+Da(ib, 1)*Hy(ib, 1);

        Eobs(t, 1) = Ez(obs, 1);
    end
    
    % fourier transform
    transform = fft(Eobs, totaltime)/totaltime;
    f=1/2/dt*linspace(0,1,totaltime/2+1)';
   
    bandstructure(k,1:range) = abs(transform(1:range,1));
    % loop to find first order peaks indicating bandgap
    % maybe use threshold to avoid zero peak
    % or use pcolor
    % figure(k);
    % plot(f(1:500,1), abs(transform(1:500, 1)));
    % pause(1);
    % for each k, find peaks in order to plot 
end

pcolor(bandstructure');
shading flat;
clim([0,0.5e-6]);
colorbar;
