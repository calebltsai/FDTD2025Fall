clear;
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
m = 2;
eta = sqrt(mu/eps);
tolerance = 1e-12;
sig_e = -log(tolerance)*(m+1)/(2*eta*pmlwidth*delta_x);

kap = 1;
alp = 0;

cw_Ez = zeros(401,1);
bw_Ez = ones(401,1);
cw_Hy = zeros(400,1);
bw_Hy = ones(400,1);
kapp_h = ones(400,1);
kapp_e = ones(401,1);

psi_Ez = zeros(401, 1);
psi_Hy = zeros(400, 1);


for x=1:pmlwidth
    iel = pmlwidth + 1 - x;
    ier = gridsize - pmlwidth + 1 + x;
    ihl = pmlwidth + 1 - x;
    ihr = gridsize - pmlwidth + x;

    sige = (sig_e*((x)/pmlwidth).^m);
    sigh = (sig_e*((x - 0.5)/pmlwidth).^m);
    
    kappe = (1+(kap-1)*(((x)/pmlwidth).^m));
    kapph = (1+(kap-1)*(((x - 0.5)/pmlwidth).^m));
    alphe = alp*(((pmlwidth - x + 1)/pmlwidth).^m);
    alphh = alp*(((pmlwidth - x + 0.5)/pmlwidth).^m);

    kapp_e(iel, 1) = kappe;
    kapp_h(ihl, 1) = kapph;
    kapp_e(ier, 1) = kappe;
    kapp_h(ihr, 1) = kapph;

    bw_Ez(iel, 1) = exp(-1*delta_t/eps*(sige/kappe+alphe));
    bw_Ez(ier, 1) = bw_Ez(iel, 1);
    cw_Ez(iel, 1) = sige*(bw_Ez(iel, 1) - 1)/(sige + kappe*alphe)/kappe;
    cw_Ez(ier, 1) = cw_Ez(iel, 1);

    bw_Hy(ihl, 1) = exp(-1*delta_t/eps*(sigh/kapph+alphh));
    bw_Hy(ihr, 1) = bw_Hy(ihl, 1);
    cw_Hy(ihl, 1) = sigh*(bw_Hy(ihl, 1) - 1)/(sigh + kapph*alphh)/kapph;
    cw_Hy(ihr, 1) = cw_Hy(ihl, 1);
end

% hold on;
% plot(linspace(1, 801, 401), bw_Ez(:,1));
% plot(linspace(2, 800, 400), bw_Hy(:,1));
% pause(5);
% hold off;


plot(ez_grid);
xlim([1 gridsize+1]);
getframe;

for t=1:pulsetime*8
    
    %update e with previous e data and 1/2 frame data of h
    for x=2:gridsize
        HyDx = (hy_grid(x,1)-hy_grid(x-1,1))/( kapp_e(x,1)*delta_x );

        psi_Ez(x,1) = bw_Ez(x,1)*psi_Ez(x,1) + cw_Ez(x,1)*HyDx;
        
        ez_grid(x,1) =  ez_grid(x,1) + delta_t/eps*( HyDx + psi_Ez(x,1));
    end
    ez_grid(gridsize/2+1,1) = ez_grid(gridsize/2+1,1) + exp(-((100-t))^2/gridsize/2)*sin(w*(t-100)*delta_t);


    %update h with previous h data and 1/2 frame data of e
    for x=1:gridsize
        EzDx = (ez_grid(x+1,1)-ez_grid(x,1))/( kapp_h(x,1)*delta_x );
        
        psi_Hy(x,1) = bw_Hy(x,1)*psi_Hy(x,1) + cw_Hy(x,1)*EzDx;
        
        hy_grid(x,1) =  hy_grid(x,1) + delta_t/mu*( EzDx + psi_Hy(x,1));
    end

    if mod(t,plotspeed) == 0
        plot(ez_grid);
        % ylim([-1 1]);
        xlim([0 gridsize + 2]);
        getframe;
        pause(0.01);
    end
end