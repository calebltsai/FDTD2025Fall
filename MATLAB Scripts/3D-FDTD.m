clear;
% % constants
w = 2*pi*1e9;
c = 299792458;
lambda = 2*pi*c/w; % recording to fit grid steps within good resolution of wavelength
mu = 4*pi*1e-7;
eps = 1/(c^2*mu);
eta = sqrt(mu/eps);

courantfactor = 0.5;
delta_x = lambda/50;
delta_y = delta_x;
delta_z = delta_x;
delta_t = courantfactor*delta_x/c;

% % grid and boundary conditions
gridsize = 400;
plotspeed = 10; % how many time evolution steps before graphed again
tolerance = 1e-18;
pmlwidth = 20;
totaltime = 200;

% % optimization perameters
m = 2.5; % gradient of easing into cpml
alp = 0; % "slows down" incident wave coming into boundary
kap = 1; % tradeoff evinescent field at boundary vs propogating field

% % initialize grid with blank frame
Ex = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
Ey = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
Ez = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
Hx = zeros(gridsize, gridsize, gridsize);
Hy = zeros(gridsize, gridsize, gridsize);
Hz = zeros(gridsize, gridsize, gridsize);

% % test recovery of snells law and transmission/reflection coefficients
% % and maybe fourier analysis of frequencies involved using
% % dimagnetic/dielectric properties
epss = ones(gridsize + 1, gridsize + 1, gridsize + 1)*eps;
mus = ones(gridsize, gridsize, gridsize)*mu;

% % setting dielectrics
% for x=1:gridsize+1
%     for y=1:gridsize
%         if (x-150).^2+(y-201).^2 <= 20
%             epss(x,y) = epss(x,y)*2;
%         end
%     end
% end
% epss(195:195+10,125:gridsize + 1-125) = 5*eps;

% % setting dimagnetics
% for x=1:gridsize+1
%     for y=1:gridsize
%     end
% end
% mus(195:195+10,125:gridsize + 1-125) = 5*mu;


% % building array of coeffs/boundary conditions
sig_max = -log(tolerance)*(m+1)/(2*eta*pmlwidth*delta_x);

% % setting E CPML consts
c_Ezx = zeros(gridsize + 1, 1);
b_Ezx = ones(gridsize + 1, 1);
c_Ezy = zeros(gridsize + 1, 1);
b_Ezy = ones(gridsize + 1, 1);
c_Eyx = zeros(gridsize + 1, 1);
b_Eyx = ones(gridsize + 1, 1);
c_Eyz = zeros(gridsize + 1, 1);
b_Eyz = ones(gridsize + 1, 1);
c_Exz = zeros(gridsize + 1, 1);
b_Exz = ones(gridsize + 1, 1);
c_Exy = zeros(gridsize + 1, 1);
b_Exy = ones(gridsize + 1, 1);

% % setting H CPML consts
c_Hyx = zeros(gridsize, 1);
b_Hyx = ones(gridsize, 1);
c_Hyz = zeros(gridsize, 1);
b_Hyz = ones(gridsize, 1);
c_Hxy = zeros(gridsize, 1);
b_Hxy = ones(gridsize, 1);
c_Hxz = zeros(gridsize, 1);
b_Hxz = ones(gridsize, 1);
c_Hzy = zeros(gridsize, 1);
b_Hzy = ones(gridsize, 1);
c_Hzx = zeros(gridsize, 1);
b_Hzx = ones(gridsize, 1);


kapp_Hyx = ones(gridsize, 1);
kapp_Hyz = ones(gridsize, 1);
kapp_Hxy = ones(gridsize, 1);
kapp_Hxz = ones(gridsize, 1);
kapp_Hzy = ones(gridsize, 1);
kapp_Hzx = ones(gridsize, 1);
kapp_Ezx = ones(gridsize + 1, 1);
kapp_Ezy = ones(gridsize + 1, 1);
kapp_Eyz = ones(gridsize + 1, 1);
kapp_Eyx = ones(gridsize + 1, 1);
kapp_Exz = ones(gridsize + 1, 1);
kapp_Exy = ones(gridsize + 1, 1);

% % arrays for pml boundary resable in square grid and uniform boundary
% % widths, coeffs increase corresponding to index
for x=1:pmlwidth
    % % initialize used indices on boundary according to gradient
    iel = pmlwidth + 1 - x;
    ier = gridsize - pmlwidth + 1 + x;
    ihl = pmlwidth + 1 - x;
    ihr = gridsize - pmlwidth + x;
    
    % % compute values used in 1:pmlwidth gradient
    % % offset by 1/2 point and non-zero initial value
    sige = (sig_max*((x)/pmlwidth).^m);
    sigh = (sig_max*((x - 0.5)/pmlwidth).^m);

    kappe = (1+(kap-1)*(((x)/pmlwidth).^m));
    kapph = (1+(kap-1)*(((x - 0.5)/pmlwidth).^m));

    alphe = alp*(((pmlwidth - x + 1)/pmlwidth).^m);
    alphh = alp*(((pmlwidth - x + 0.5)/pmlwidth).^m);

    % % kappa gradient used in time evolution
    kapp_Ezx(iel, 1) = kappe;
    kapp_Ezx(ier, 1) = kappe;
    kapp_Ezy(iel, 1) = kappe;
    kapp_Ezy(ier, 1) = kappe;
    kapp_Eyz(iel, 1) = kappe;
    kapp_Eyz(ier, 1) = kappe;
    kapp_Eyx(iel, 1) = kappe;
    kapp_Eyx(ier, 1) = kappe;
    kapp_Exy(iel, 1) = kappe;
    kapp_Exy(ier, 1) = kappe;
    kapp_Exz(iel, 1) = kappe;
    kapp_Exz(ier, 1) = kappe;

    kapp_Hzx(iel, 1) = kapph;
    kapp_Hzx(ier, 1) = kapph;
    kapp_Hzy(iel, 1) = kapph;
    kapp_Hzy(ier, 1) = kapph;
    kapp_Hyz(iel, 1) = kapph;
    kapp_Hyz(ier, 1) = kapph;
    kapp_Hyx(iel, 1) = kapph;
    kapp_Hyx(ier, 1) = kapph;
    kapp_Hxy(iel, 1) = kapph;
    kapp_Hxy(ier, 1) = kapph;
    kapp_Hxz(iel, 1) = kapph;
    kapp_Hxz(ier, 1) = kapph;
    
    % % CPML b and c coefficients for psi functions in Ezx, Ezy, Hy, Hx
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

% % initialize psi functions CPML
psi_Ezx = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
psi_Ezy = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
psi_Eyx = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
psi_Eyz = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
psi_Exy = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
psi_Exz = zeros(gridsize + 1, gridsize + 1, gridsize + 1);
psi_Hzx = zeros(gridsize, gridsize, gridsize);
psi_Hzy = zeros(gridsize, gridsize, gridsize);
psi_Hyx = zeros(gridsize, gridsize, gridsize);
psi_Hyz = zeros(gridsize, gridsize, gridsize);
psi_Hxy = zeros(gridsize, gridsize, gridsize);
psi_Hxz = zeros(gridsize, gridsize, gridsize);

% % capture initial frame 2D
f = figure(1);
imagesc(Ez(:, :, gridsize/2+1));
pcolor(Ez(:, :, gridsize/2+1)');
shading flat;
view(2);
daspect([1,1,1]);
colorbar;
clim([-0.2,0.2]);
% initialize animation
animation = [getframe(f)];
pause(100);

% % capture initial frame 3D
% figure(1);
% surf(ez_grid, 'edgecolor', 'none');
% daspect([1,1,1]);
% colorbar;
% clim([-0.1,0.1]);
% zlim([-0.1,0.1]);
% getframe();

% % core time evolution simulation
global time;
time = 0;

t = timer('TimerFcn', {@ProgressCheck, totaltime}, 'ExecutionMode', 'fixedRate');
t.TasksToExecute = Inf;
function ProgressCheck(obj, event, buh)
    global time;
    if(time == 0)
        disp("rendering animation:");
    elseif(time == buh)
        disp("rendering progress: 100%");
        disp("begin animation playback");
    else
        disp("rendering progress: "+num2str(100*time/buh)+"%");
    end
end
start(t);

% % time step evolution
while time < totaltime
    % % update e with previous e data and 1/2 frame data of h
    for x=2:gridsize
        for y=2:gridsize
            for z=2:gridsize
                epsilon = eps(x, y, z);

                HxDy = (Hx(x, y, z)-Hx(x, y-1, z))/(kapp_Hxy(y, 1)*delta_y);
                HxDz = (Hx(x, y, z)-Hx(x, y, z-1))/(kapp_Hxz(z, 1)*delta_z);
                
                HyDx = (Hy(x, y, z)-Hy(x-1, y, z))/(kapp_Hxy(x, 1)*delta_x);
                HyDz = (Hy(x, y, z)-Hy(x, y, z-1))/(kapp_Hxz(z, 1)*delta_z);

                HzDx = (Hz(x, y, z)-Hz(x-1, y, z))/(kapp_Hxy(x, 1)*delta_x);
                HzDy = (Hz(x, y, z)-Hz(x, y-1, z))/(kapp_Hxz(y, 1)*delta_y);
                
                psi_Ezx(x, y, z) = b_Ezx(x, 1)*psi_Ezx(x, y, z) + c_EzAx(x, 1)*HyDx;
                psi_Ezy(x, y, z) = b_Ezy(y, 1)*psi_Ezy(x, y, z) + c_Ezy(y, 1)*HxDy;
                
                Ez(x, y, z) =  Ez(x, y, z) + delta_t/epsilon*( HyDx - HxDy + psi_Ezx(x, y, z) - psi_Ezy(x, y, z) );
                Ey(x, y, z) =  Ey(x, y, z) - delta_t/epsilon*( HzDx - HxDz + psi_Eyx(x, y, z) - psi_Eyz(x, y, z) );
                Ex(x, y, z) =  Ex(x, y, z) + delta_t/epsilon*( HzDy - HyDz + psi_Exy(x, y, z) - psi_Exz(x, y, z) );
            end
        end
    end
    % % point source at center
    Ez(gridsize/2 + 1, gridsize/2 + 1, gridsize/2 + 1) = exp(-((time-200)/100)^2)*sin(w*(time-200)*delta_t) + Ez(gridsize/2, gridsize/2 + 1, gridsize/2 + 1);
    

    for x=1:gridsize
        for y=1:gridsize
            for z=1:gridsize
                % % total/scattered fields
                EzDy = (Ez(x, y+1, z)-Ez(x, y, z))/(kapp_ey(y, 1)*delta_y);
                EzDx = (Ez(x+1, y, z)-Ez(x, y, z))/(kapp_ex(x, 1)*delta_x);
    
                EyDz = (Ey(x, y, z+1)-Ey(x, y, z))/(kapp_ez(z,1)*delta_z);
                EyDx = (Ey(x+1, y, z)-Ey(x, y, z))/(kapp_ex(x,1)*delta_x);
    
                ExDy = (Ex(x, y+1, z)-Ex(x, y, z))/(kapp_ey(y,1)*delta_y);
                ExDz = (Ex(x, y, z+1)-Ex(x, y, z))/(kapp_ez(z,1)*delta_z);
    
                psi_Hxy(x, y, z) = b_Hxy(y, 1)*psi_Hxy(x,y) + c_Hxy(y, 1)*EzDy;
                psi_Hxz(x, y, z) = b_Hxy(z, 1)*psi_Hxy(x,y) + c_Hxy(y, 1)*EzDy;

                psi_Hyx(x,y) = b_Hyx(x, 1)*psi_Hyx(x,y) + c_Hyx(x, 1)*EzDx;
                psi_Hyx(x,y) = b_Hyx(x, 1)*psi_Hyx(x,y) + c_Hyx(x, 1)*EzDx;
    
                Hy(x,y) = Hy(x,y) + delta_t/mus(x,y)*( EzDx + psi_Hyx(x, y));
                Hx(x,y) = Hx(x,y) - delta_t/mus(x,y)*( EzDy + psi_Hxy(x, y));
            end
        end
    end



   
    if mod(time, plotspeed)==0
        
        % % change slice of Ez shown
        % % 2D
        imagesc(Ez);
        pcolor(Ez');
        shading flat;
        view(2);
        daspect([1,1,1]);
        colorbar;
        % clim([-0.4,0.4]);
        
 
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
pause(1);
% animation playback
stop(t);
delete(t);

% movie(animation, 1, 10);