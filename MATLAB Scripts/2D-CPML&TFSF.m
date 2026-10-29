clear;
% constants
w = 2*pi*1e9;
w_excitation = 2*pi*7.456e8; % excitation frequency derived from fourier analysis in dielectric material
c = 299792458;
lambda = 2*pi*c/w; % recording to fit grid steps within good resolution of wavelength
mu = 4*pi*1e-7;
eps = 1/(c^2*mu);
eta = sqrt(mu/eps);

courantfactor = 0.5;
delta_x = lambda/50;
delta_y = delta_x;
delta_t = courantfactor*delta_x/c;

% grid and boundary conditions
gridsize = 400;
plotspeed = 10; % how many time evolution steps before graphed again
tolerance = 1e-18;
pmlwidth = 20;
totaltime = 3000;
time = 0;

ie_sfleft = 100;
ie_sfright = 301;

% optimization perameters
m = 2.5; % gradient of easing into cpml
alp = 0; % "slows down" incident wave coming into boundary
kap = 1; % tradeoff evinescent field at boundary vs propogating field

% initialize grid for TMz mode (with +1 in egrid) with blank frame
% e_obs = zeros(gridsize + 1, 1);
inc_ez = zeros(gridsize + 1,1);
inc_hy = zeros(gridsize,1);

ez_grid = zeros(gridsize + 1);
hy_grid = zeros(gridsize);
hx_grid = zeros(gridsize);

% test recovery of snells law and transmission/reflection coefficients
% and maybe fourier analysis of frequencies involved using
% dimagnetic/dielectric properties
epss = ones(gridsize + 1)*eps;
mus = ones(gridsize)*mu;

% % setting dielectrics
% for x=1:gridsize+1
%     for y=1:gridsize
%         if (x-150).^2+(y-201).^2 <= 20
%             epss(x,y) = epss(x,y)*2;
%         end
%     end
% end
epss(197:197+5,150:gridsize + 1-150) = 5*eps;

% % setting dimagnetics
% for x=1:gridsize+1
%     for y=1:gridsize
%     end
% end
% mus(195:195+10,125:gridsize + 1-125) = 5*mu;


% % for testing numerical error and optimization through reflective error of pml
% large_ez = zeros(gridsize*3+1);
% large_hy = zeros(gridsize*3);
% large_hx = zeros(gridsize*3);

% building array of coeffs/boundary conditions
sig_max = -log(tolerance)*(m+1)/(2*eta*pmlwidth*delta_x);

c_Ezx = zeros(gridsize + 1, 1); % generally would use 2D grid to encompass all material
c_Ezy = zeros(gridsize + 1, 1);
b_Ezx = ones(gridsize + 1, 1);
b_Ezy = ones(gridsize + 1, 1);

c_Hyx = zeros(gridsize, 1);
c_Hxy = zeros(gridsize, 1);
b_Hyx = ones(gridsize, 1);
b_Hxy = ones(gridsize, 1);

kapp_hy = ones(gridsize, 1);
kapp_hx = ones(gridsize, 1);
kapp_ex = ones(gridsize + 1, 1);
kapp_ey = ones(gridsize + 1, 1);

% arrays for pml boundary resable in square grid and uniform boundary
% widths, coeffs increase corresponding to index
for x=1:pmlwidth
    % initialize used indices on boundary according to gradient
    iel = pmlwidth + 1 - x;
    ier = gridsize - pmlwidth + 1 + x;
    ihl = pmlwidth + 1 - x;
    ihr = gridsize - pmlwidth + x;
    
    % compute values used in 1:pmlwidth gradient
    % offset by 1/2 point and non-zero initial value
    sige = (sig_max*((x)/pmlwidth).^m);
    sigh = (sig_max*((x - 0.5)/pmlwidth).^m);

    kappe = (1+(kap-1)*(((x)/pmlwidth).^m));
    kapph = (1+(kap-1)*(((x - 0.5)/pmlwidth).^m));

    alphe = alp*(((pmlwidth - x + 1)/pmlwidth).^m);
    alphh = alp*(((pmlwidth - x + 0.5)/pmlwidth).^m);

    % kappa gradient used in time evolution
    kapp_ex(iel, 1) = kappe;
    kapp_ex(ier, 1) = kappe;
    kapp_ey(iel, 1) = kappe;
    kapp_ey(ier, 1) = kappe;

    kapp_hx(iel, 1) = kapph;
    kapp_hx(ier, 1) = kapph;
    kapp_hy(iel, 1) = kapph;
    kapp_hy(ier, 1) = kapph;
    
    % CPML b and c coefficients for psi functions in Ezx, Ezy, Hy, Hx
    b_Ezx(iel, 1) = exp(-1*delta_t/eps*(sige/kappe+alphe));
    b_Ezx(ier, 1) = b_Ezx(iel, 1);
    c_Ezx(iel, 1) = sige*(b_Ezx(iel, 1) - 1)/(sige + kappe*alphe)/kappe;
    c_Ezx(ier, 1) = c_Ezx(iel, 1);

    b_Ezy(iel, 1) = exp(-1*delta_t/eps*(sige/kappe+alphe));
    b_Ezy(ier, 1) = b_Ezy(iel, 1);
    c_Ezy(iel, 1) = sige*(b_Ezy(iel, 1) - 1)/(sige + kappe*alphe)/kappe;
    c_Ezy(ier, 1) = c_Ezy(iel, 1);

    b_Hyx(ihl, 1) = exp(-1*delta_t/eps*(sigh/kapph+alphh));
    b_Hyx(ihr, 1) = b_Hyx(ihl, 1);
    c_Hyx(ihl, 1) = sigh*(b_Hyx(ihl, 1) - 1)/(sigh + kapph*alphh)/kapph;
    c_Hyx(ihr, 1) = c_Hyx(ihl, 1);
    b_Hxy(ihl, 1) = exp(-1*delta_t/eps*(sigh/kapph+alphh));
    b_Hxy(ihr, 1) = b_Hxy(ihl, 1);
    c_Hxy(ihl, 1) = sigh*(b_Hxy(ihl, 1) - 1)/(sigh + kapph*alphh)/kapph;
    c_Hxy(ihr, 1) = c_Hxy(ihl, 1);

end

% % check geometry/gradient of coefficient
% hold on;
% plot(kapp_ex);
% pause(5);
% plot(kapp_ey);
% pause(5);
% plot(kapp_hy);
% pause(5);
% plot(kapp_hx);
% pause(20);

% hold on;
% plot(linspace(1, 801, gridsize + 1), c_Ezx(:,1));
% plot(linspace(2, 800, gridsize), c_Hyx(:,1));
% plot(linspace(2, 800, gridsize), c_Hxy(:,1));
% pause(50);
% hold off;

psi_Ezx = zeros(gridsize + 1);
psi_Ezy = zeros(gridsize + 1);
psi_Hxy = zeros(gridsize);
psi_Hyx = zeros(gridsize);

psi_Ezx1 = zeros(gridsize + 1,1);
psi_Hyx1 = zeros(gridsize,1);

% fourier analysis variables
nfreq = 2000;
w_start = 2*pi*1e8;
w_end = 2*pi*3e9;
ez_dft = zeros(nfreq, 1);
omega = linspace(w_start, w_end, nfreq)';
ez_grid_dft = zeros(gridsize+1);

% capture initial frame 2D
f = figure(1);
imagesc(ez_grid);
pcolor(ez_grid');
shading flat;
view(2);
daspect([1,1,1]);
colorbar;
clim([-0.2,0.2]);
% initialize animation
animation = [getframe(f)];

% % capture initial frame 3D
% figure(1);
% surf(ez_grid, 'edgecolor', 'none');
% daspect([1,1,1]);
% colorbar;
% clim([-0.1,0.1]);
% zlim([-0.1,0.1]);
% getframe();

% core time evolution simulation

% time step evolution
while time < totaltime
    % update e with previous e data and 1/2 frame data of h
    for x=2:gridsize
        HyDx1 = (inc_hy(x,1)-inc_hy(x-1,1))/(kapp_hx(x,1)*delta_x);
        psi_Ezx1(x,1) = b_Ezx(x, 1)*psi_Ezx1(x,1) + c_Ezx(x, 1)*HyDx1;
        inc_ez(x,1) =  inc_ez(x,1) + delta_t/eps*( HyDx1 + psi_Ezx1(x,1) );
        for y=2:gridsize
            HxDy = (hx_grid(x,y)-hx_grid(x,y-1))/(kapp_hy(y,1)*delta_y);
            HyDx = (hy_grid(x,y)-hy_grid(x-1,y))/(kapp_hx(x,1)*delta_x);

            psi_Ezx(x,y) = b_Ezx(x, 1)*psi_Ezx(x,y) + c_Ezx(x, 1)*HyDx;
            psi_Ezy(x,y) = b_Ezy(y, 1)*psi_Ezy(x,y) + c_Ezy(y, 1)*HxDy;
            
            ez_grid(x,y) =  ez_grid(x,y) + delta_t/epss(x,y)*( HyDx - HxDy + psi_Ezx(x,y) - psi_Ezy(x,y) );
        end
    end
    % soft source
    % ez_grid(gridsize/2,gridsize/2) = exp(-((time-200)/100)^2)*sin(w*(time-200)*delta_t) + ez_grid(gridsize/2,gridsize/2);
    % inc_ez(gridsize/2,gridsize/2) = exp(-((time-200)/100)^2)*sin(w*(time-200)*delta_t) + inc_ez(gridsize/2,gridsize/2);

    % scattered boundary correction
    for i=ie_sfleft:ie_sfright
        ez_grid(ie_sfleft,i) = ez_grid(ie_sfleft,i) - delta_t/(delta_x*epss(ie_sfleft,i))*inc_hy(ie_sfleft-1, 1);
        ez_grid(ie_sfright,i) = ez_grid(ie_sfright,i) + delta_t/(delta_x*epss(ie_sfright,i))*inc_hy(ie_sfright, 1);
    end
    inc_ez(40,1) = exp(-((time-100)/50)^2)*sin(w_excitation*(time-100)*delta_t) + inc_ez(40,1);

    % ez_grid(100:300,gridsize/2-15) = 0;
    % ez_grid(100:300,gridsize/2+15) = 0;

    % e_obs(t, 1) = ez_grid(gridsize/2 + 20, gridsize/2);

    % update h with previous h data and 1/2 frame data of e
    for x=1:gridsize
        EzDx1 = (inc_ez(x+1,1)-inc_ez(x,1))/(kapp_ex(x,1)*delta_x);
        psi_Hyx1(x,1) = b_Hyx(x, 1)*psi_Hyx1(x,1) + c_Hyx(x, 1)*EzDx1;
        inc_hy(x,1) = inc_hy(x,1) + delta_t/mu*( EzDx1 + psi_Hyx1(x, 1) );
        for y=1:gridsize
            % total/scattered fields
            EzDy = (ez_grid(x,y+1)-ez_grid(x,y))/(kapp_ey(y,1)*delta_y);
            EzDx = (ez_grid(x+1,y)-ez_grid(x,y))/(kapp_ex(x,1)*delta_x);

            psi_Hxy(x,y) = b_Hxy(y, 1)*psi_Hxy(x,y) + c_Hxy(y, 1)*EzDy;
            psi_Hyx(x,y) = b_Hyx(x, 1)*psi_Hyx(x,y) + c_Hyx(x, 1)*EzDx;

            hy_grid(x,y) = hy_grid(x,y) + delta_t/mus(x,y)*( EzDx + psi_Hyx(x, y));
            hx_grid(x,y) = hx_grid(x,y) - delta_t/mus(x,y)*( EzDy + psi_Hxy(x, y));
        end
    end

    % scattered boundary correction
    for j=ie_sfleft:ie_sfright
        hy_grid(ie_sfright,j) = hy_grid(ie_sfright,j) + delta_t/(delta_x*mus(ie_sfright,j))*(inc_ez(ie_sfright,1));
        hy_grid(ie_sfleft-1,j) = hy_grid(ie_sfleft-1,j) - delta_t/(delta_x*mus(ie_sfleft-1,j))*(inc_ez(ie_sfleft,1));
        hx_grid(j,ie_sfright) = hx_grid(j,ie_sfright) - delta_t/(delta_x*mus(j,ie_sfright))*(inc_ez(j,1));
        hx_grid(j,ie_sfleft-1) = hx_grid(j,ie_sfleft-1) + delta_t/(delta_x*mus(j,ie_sfleft-1))*(inc_ez(j,1));
    end
    
    % fourier analysis of wavelength
    for iff=1:nfreq
        ez_dft(iff,1) = ez_dft(iff,1) + ez_grid(200, 200)*exp(1j*omega(iff,1)*time*delta_t);
    end

    ez_grid_dft(2:gridsize,2:gridsize) = ez_grid_dft(2:gridsize,2:gridsize) + ez_grid(2:gridsize,2:gridsize)*exp(1j*w_excitation*time*delta_t);
    
   
    if mod(time, plotspeed)==0
        % % show outline of material
        % ez_grid(195, 125:gridsize+1-125) = ez_grid(195, 125:gridsize+1-125) + 1;
        % ez_grid(195+10, 125:gridsize+1-125) = ez_grid(195+10, 125:gridsize+1-125) + 1;
        % ez_grid(195+1:195+10-1, 125) = ez_grid(195+1:195+10-1, 125) + 1;
        % ez_grid(195+1:195+10-1,gridsize + 1 - 125) = ez_grid(195+1:195+10-1, gridsize + 1 - 125) + 1;

        % 2D
        imagesc(ez_grid);
        pcolor(ez_grid');
        shading flat;
        view(2);
        daspect([1,1,1]);
        colorbar;
        % clim([-0.01,0.01]);
        
        % % reset e grid after showing boundary
        % ez_grid(195, 125:gridsize + 1 - 125) = ez_grid(195, 125:gridsize + 1 - 125) - 1;
        % ez_grid(195+10, 125:gridsize + 1 - 125) = ez_grid(195+10, 125:gridsize + 1 - 125) - 1;
        % ez_grid(195+1:195+10-1, 125) = ez_grid(195+1:195+10-1, 125) - 1;
        % ez_grid(195+1:195+10-1, gridsize + 1 - 125) = ez_grid(195+1:195+10-1, gridsize + 1 - 125) - 1;
 
        % % 3D
        % figure(1);
        % set(gcf,'Visible', 'off');
        % imagesc(ez_grid);
        % surf(ez_grid, 'edgecolor', 'none');
        % daspect([1,1,1]);
        % zlim([-0.1 0.1]);
        % colorbar;
        % clim([-0.1,0.1]);

        animation(time/plotspeed + 1, 1) = getframe();

    end
    time = time + 1;
end

figure(2);
plot(omega/2/pi, abs(ez_dft));

figure(3);
pcolor(abs(ez_grid_dft)');
shading flat;
daspect([1 1 1]);
colorbar;

% animation playback

% movie(animation, 1, 10);