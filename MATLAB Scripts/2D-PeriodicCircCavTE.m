clear;
% constants
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
dy = dx;

% grid and boundary conditions
ib = 32;
jb = 32;
totaltime = (2^10)*(2^4);


Ez = zeros(ib+1);
Hy = zeros(ib+1);
Hx = zeros(ib+1);

epss = ones(ib+1)*epsilon;
mus = ones(ib+1)*mu;

% setting dielectrics
for x=1:ib
    for y=1:ib
        if (x-ib/2).^2+(y-ib/2).^2 <= (ib/4).^2
            epss(x,y) = epss(x,y)*(2.^2);
        end
    end
end

% % setting dimagnetics
% for x=1:gridsize+1
%     for y=1:gridsize
%     end
% end


% handle the array of k vectors in reduced brillouin zone used for bloch condition
numK = 10;
kxs = [linspace(0, pi/(dx*ib), numK+1), ones(1, numK+1)*pi/(dx*ib), linspace(pi/(dy*ib), 0, numK+1)]' ;
kys = [zeros(1,numK+1), linspace(0, pi/(dy*ib), numK+1), linspace(pi/(dy*ib), 0, numK+1)]';


range = 500;

bandstructure = zeros(numK, range);

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
    % f=1/2/dt*linspace(0,1,totaltime/2+1)';

    bandstructure(k,1:range) = abs(transform(1:range,1));
end

figure(1);
pcolor(bandstructure');
shading flat;
clim([0,0.5e-6]);
colorbar;