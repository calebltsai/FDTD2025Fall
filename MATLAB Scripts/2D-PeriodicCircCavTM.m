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
ib = 200;
jb = 200;
totaltime = (2^10)*(2^4);


Hz = zeros(ib+1);
Ey = zeros(ib+1);
Ex = zeros(ib+1);

epss = ones(ib+1)*epsilon;
mus = ones(ib+1)*mu;

% setting dielectrics
for x=1:ib
    for y=1:ib
        if (x-ib/2).^2+(y-ib/2).^2 >= (ib/3).^2
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
    h_obs = zeros(totaltime, 1);
    % time step evolution
    for time=1:totaltime
        % update e with previous e data and 1/2 frame data of h
        for x=2:ib+1
            for y=2:jb+1
                HzDy = ( Hz(x,y)-Hz(x,y-1) )/dy;
                HzDx = ( Hz(x,y)-Hz(x-1,y) )/dx;
    
                Ex(x,y) = Ex(x,y) + dt/epss(x,y)*( HzDy );
                Ey(x,y) = Ey(x,y) - dt/epss(x,y)*( HzDx );
            end
        end
        Ey(1,:)=Ey(ib+1,:)*exp(1i*kxs(k, 1)*ib*dx);
        Ex(:,1)=Ex(:,jb+1)*exp(1j*kys(k, 1)*jb*dy);

        % scattered boundary correction
        % update h with previous h data and 1/2 frame data of e

        for x=1:ib
            for y=1:jb
                % total/scattered fields
                ExDy = ( Ex(x,y+1)-Ex(x,y) )/dy;
                EyDx = ( Ey(x+1,y)-Ey(x,y) )/dx;
    
                Hz(x,y) = Hz(x,y) + dt/mus(x,y)*(  ExDy - EyDx );
            end
        end
        Hz(:,jb+1) = Hz(:,1)*exp(-1j*kys(k,1)*jb*dy);
        Hz(ib+1,:) = Hz(1,:)*exp(-1i*kxs(k,1)*ib*dx);

        Hz(ib/4, jb/4) = exp(-(time-10)^2/((2)^2/10))*sin(omega*dt*(time-10)) + Hz(ib/4,jb/4);

        h_obs(time, 1) = Hz(ib/2, jb/2);
        pcolor(Hz);
        shading flat;
        colorbar;
        pause(0.01);
    end

    transform = fft(h_obs, totaltime)/totaltime;
    % f=1/2/dt*linspace(0,1,totaltime/2+1)';

    bandstructure(k,1:range) = abs(transform(1:range,1));
end

figure(2);
pcolor(bandstructure');
shading flat;
clim([0,0.2e-6]);
colorbar;