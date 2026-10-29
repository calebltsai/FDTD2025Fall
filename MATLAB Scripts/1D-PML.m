clear();
w = 2*pi*1e9;
c = 299792458;
lambda = 2*pi*c/w;
mu = 4*pi*1e-7;
eps = 1/(c^2*mu);
delta_x = lambda/40;
delta_t = delta_x/c*0.5;

% convenient params to change
gridsize = 400;
pulsetime = 100;
plotspeed = 4;


%initialize with blank frame
ez_grid = zeros(gridsize+1,1);
hy_grid = zeros(gridsize,1);

pmlwidth = 20;
m = 2.5;

eta = sqrt(mu/eps);
sig_e = -log(1e-20)*(m+1)/(2*eta*pmlwidth*delta_x);
sig_m = sig_e*eta*eta;
% in pml sig_e*((gridsize-x)/pmlwidth).^3

plot(ez_grid);
xlim([0 gridsize]);
getframe;

ca = ones(401,1)*delta_t/eps;
cb = ones(401,1)*eps/delta_t;
% same for h -> da, db

for x=2:pmlwidth+1
    ca(x,1) = 1/(eps/delta_t+(sig_e*((pmlwidth+1-x)/pmlwidth).^m)/2);
    cb(x,1) = eps/delta_t-(sig_e*((pmlwidth+1-x)/pmlwidth).^m)/2;
end
for x=gridsize-pmlwidth+1:gridsize
    ca(x,1) = 1/(eps/delta_t+(sig_e*((x-gridsize-1+pmlwidth)/pmlwidth).^m)/2);
    cb(x,1) = (eps/delta_t-(sig_e*((x-gridsize-1+pmlwidth)/pmlwidth).^m)/2);
end

for t=1:pulsetime*8
    
    %update e with previous e data and 1/2 frame data of h
    % ez_grid(1,1) = sin(w*t*delta_t);
    for x=2:gridsize
        ez_grid(x,1) = ca(x,1)*(cb(x,1)*ez_grid(x,1)+(hy_grid(x,1)-hy_grid(x-1,1))/delta_x);
    end
    ez_grid(gridsize/2+1,1) = ez_grid(gridsize/2+1,1) + exp(-(100-t)^2/gridsize/2)*sin(w*(t-100)*delta_t);

    

    %update h with previous h data and 1/2 frame data of e
    for x=1:pmlwidth
        hy_grid(x,1) = 1/(mu/delta_t+sig_m*((pmlwidth-x+0.5)/pmlwidth).^m/2)*(((mu/delta_t-sig_m*((pmlwidth-x+0.5)/pmlwidth).^m/2))*hy_grid(x,1)+(ez_grid(x+1,1)-ez_grid(x,1))/delta_x);
    end
    for x=pmlwidth+1:gridsize-pmlwidth
        hy_grid(x,1) = delta_t/(delta_x*mu)*(ez_grid(x+1,1)-ez_grid(x,1))+hy_grid(x,1);      
    end
    for x=gridsize-pmlwidth+1:gridsize
        hy_grid(x,1) = 1/(mu/delta_t+sig_m*((x-gridsize+pmlwidth-0.5)/pmlwidth).^m/2)*(((mu/delta_t-sig_m*((x+pmlwidth-gridsize-0.5)/pmlwidth).^m/2))*hy_grid(x,1)+(ez_grid(x+1,1)-ez_grid(x,1))/delta_x);
    end

    if mod(t,plotspeed) == 0
        plot(ez_grid);
        xlim([0 gridsize]);
        getframe;
    end
end