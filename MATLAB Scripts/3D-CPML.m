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
gridsize = 100;
plotspeed = 4; % how many time evolution steps before graphed again
tolerance = 1e-18;
pmlwidth = 20;
totaltime = 800;

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
% epss(gridsize/2-5:gridsize/2+5,gridsize/2-5:gridsize/2+5,125:gridsize + 1-125) = 5*eps;

% % setting dimagnetics
% for x=1:gridsize+1
%     for y=1:gridsize
%     end
% end
% mus(gridsize/2-5:gridsize/2+5,gridsize/2-5:gridsize/2+5,125:gridsize + 1-125, ) = 5*mu;


% % building array of coeffs/boundary conditions
sig_max = -log(tolerance)*(m+1)/(2*eta*pmlwidth*delta_x);

% % setting E CPML consts
c_Ex = zeros(gridsize + 1, 1);
b_Ex = ones(gridsize + 1, 1);
c_Ey = zeros(gridsize + 1, 1);
b_Ey = ones(gridsize + 1, 1);
c_Ez = zeros(gridsize + 1, 1);
b_Ez = ones(gridsize + 1, 1);

% % setting H CPML consts
c_Hx = zeros(gridsize, 1);
b_Hx = ones(gridsize, 1);
c_Hy = zeros(gridsize, 1);
b_Hy = ones(gridsize, 1);
c_Hz = zeros(gridsize, 1);
b_Hz = ones(gridsize, 1);


kapp_Hx = ones(gridsize, 1);
kapp_Hy = ones(gridsize, 1);
kapp_Hz = ones(gridsize, 1);
kapp_Ex = ones(gridsize + 1, 1);
kapp_Ey = ones(gridsize + 1, 1);
kapp_Ez = ones(gridsize + 1, 1);

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
    kapp_Ex(iel, 1) = kappe;
    kapp_Ex(ier, 1) = kappe;
    kapp_Ey(iel, 1) = kappe;
    kapp_Ey(ier, 1) = kappe;
    kapp_Ez(iel, 1) = kappe;
    kapp_Ez(ier, 1) = kappe;

    kapp_Hx(ihl, 1) = kapph;
    kapp_Hx(ihr, 1) = kapph;
    kapp_Hy(ihl, 1) = kapph;
    kapp_Hy(ihr, 1) = kapph;
    kapp_Hz(ihl, 1) = kapph;
    kapp_Hz(ihr, 1) = kapph;
    
    % % CPML b and c coefficients for psi functions
    b_Ex(iel, 1) = exp(-1*delta_t/eps*(sige/kappe+alphe));
    b_Ex(ier, 1) = b_Ex(iel, 1);
    c_Ex(iel, 1) = sige*(b_Ex(iel, 1) - 1)/(sige + kappe*alphe)/kappe;
    c_Ex(ier, 1) = c_Ex(iel, 1);

    b_Ey(iel, 1) = exp(-1*delta_t/eps*(sige/kappe+alphe));
    b_Ey(ier, 1) = b_Ey(iel, 1);
    c_Ey(iel, 1) = sige*(b_Ey(iel, 1) - 1)/(sige + kappe*alphe)/kappe;
    c_Ey(ier, 1) = c_Ey(iel, 1);

    b_Ez(iel, 1) = exp(-1*delta_t/eps*(sige/kappe+alphe));
    b_Ez(ier, 1) = b_Ez(iel, 1);
    c_Ez(iel, 1) = sige*(b_Ez(iel, 1) - 1)/(sige + kappe*alphe)/kappe;
    c_Ez(ier, 1) = c_Ez(iel, 1);

    b_Hx(ihl, 1) = exp(-1*delta_t/eps*(sigh/kapph+alphh));
    b_Hx(ihr, 1) = b_Hx(ihl, 1);
    c_Hx(ihl, 1) = sigh*(b_Hx(ihl, 1) - 1)/(sigh + kapph*alphh)/kapph;
    c_Hx(ihr, 1) = c_Hx(ihl, 1);

    b_Hy(ihl, 1) = exp(-1*delta_t/eps*(sigh/kapph+alphh));
    b_Hy(ihr, 1) = b_Hy(ihl, 1);
    c_Hy(ihl, 1) = sigh*(b_Hy(ihl, 1) - 1)/(sigh + kapph*alphh)/kapph;
    c_Hy(ihr, 1) = c_Hy(ihl, 1);

    b_Hz(ihl, 1) = exp(-1*delta_t/eps*(sigh/kapph+alphh));
    b_Hz(ihr, 1) = b_Hz(ihl, 1);
    c_Hz(ihl, 1) = sigh*(b_Hz(ihl, 1) - 1)/(sigh + kapph*alphh)/kapph;
    c_Hz(ihr, 1) = c_Hz(ihl, 1);
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

% % capture initial frame 2D slice
f = figure(1);
slice(Ez,gridsize/4,[],[]);
hold on; 
% slice(Ez,gridsize + 1 - gridsize/4,[],[]);
% slice(Ez,[],gridsize/4,[]);
% slice(Ez,[],gridsize + 1 - gridsize/4,[]);
% slice(Ez,[],[],gridsize/4);
% slice(Ez,[],[],gridsize + 1 - gridsize/4);
hold off;
shading flat;
colorbar;
daspect([1 1 1]);
clim([-0.005 0.005]);
% initialize animation
animation = [getframe(f)];
disp("rendering animation:");

% % core time evolution simulation
global time;
time = 0;

t = timer('TimerFcn', {@ProgressCheck, totaltime}, 'ExecutionMode', 'fixedRate');
t.TasksToExecute = Inf;
t.period = 5;
function ProgressCheck(~, ~, buh)
    global time;
    if(time >= buh)
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
                epsilon = epss(x, y, z);
                % update CPML expansion
                HxDy = (Hx(x, y, z)-Hx(x, y-1, z))/(kapp_Hy(y, 1)*delta_y);
                HxDz = (Hx(x, y, z)-Hx(x, y, z-1))/(kapp_Hz(z, 1)*delta_z);
                
                HyDx = (Hy(x, y, z)-Hy(x-1, y, z))/(kapp_Hx(x, 1)*delta_x);
                HyDz = (Hy(x, y, z)-Hy(x, y, z-1))/(kapp_Hz(z, 1)*delta_z);

                HzDx = (Hz(x, y, z)-Hz(x-1, y, z))/(kapp_Hx(x, 1)*delta_x);
                HzDy = (Hz(x, y, z)-Hz(x, y-1, z))/(kapp_Hy(y, 1)*delta_y);

                psi_Exy(x, y, z) = b_Ey(x, 1)*psi_Exy(x, y, z) + c_Ey(x, 1)*HzDy;
                psi_Exz(x, y, z) = b_Ez(y, 1)*psi_Exz(x, y, z) + c_Ez(y, 1)*HyDz;

                psi_Eyx(x, y, z) = b_Ex(x, 1)*psi_Eyx(x, y, z) + c_Ex(x, 1)*HzDx;
                psi_Eyz(x, y, z) = b_Ez(z, 1)*psi_Eyz(x, y, z) + c_Ez(z, 1)*HxDz;

                psi_Ezx(x, y, z) = b_Ex(x, 1)*psi_Ezx(x, y, z) + c_Ex(x, 1)*HyDx;
                psi_Ezy(x, y, z) = b_Ey(y, 1)*psi_Ezy(x, y, z) + c_Ey(y, 1)*HxDy;
                
                Ez(x, y, z) =  Ez(x, y, z) + delta_t/epsilon*( HyDx - HxDy + psi_Ezx(x, y, z) - psi_Ezy(x, y, z) );
                Ey(x, y, z) =  Ey(x, y, z) - delta_t/epsilon*( HzDx - HxDz + psi_Eyx(x, y, z) - psi_Eyz(x, y, z) );
                Ex(x, y, z) =  Ex(x, y, z) + delta_t/epsilon*( HzDy - HyDz + psi_Exy(x, y, z) - psi_Exz(x, y, z) );
            end
        end
    end
    % % point source at center
    Ez(gridsize/2 + 1, gridsize/2 + 1, gridsize/2 + 1) = 25*exp(-((time-(totaltime/5))/(totaltime/5))^2)*sin(w*(time-(totaltime/5))*delta_t) + Ez(gridsize/2, gridsize/2 + 1, gridsize/2 + 1);
    
    for x=1:gridsize
        for y=1:gridsize
            for z=1:gridsize
                tempmu = mus(x, y, z);
                % update CPML expansion
                EzDy = (Ez(x, y+1, z)-Ez(x, y, z))/(kapp_Ey(y, 1)*delta_y);
                EzDx = (Ez(x+1, y, z)-Ez(x, y, z))/(kapp_Ex(x, 1)*delta_x);
    
                EyDz = (Ey(x, y, z+1)-Ey(x, y, z))/(kapp_Ez(z, 1)*delta_z);
                EyDx = (Ey(x+1, y, z)-Ey(x, y, z))/(kapp_Ex(x, 1)*delta_x);
    
                ExDy = (Ex(x, y+1, z)-Ex(x, y, z))/(kapp_Ey(y, 1)*delta_y);
                ExDz = (Ex(x, y, z+1)-Ex(x, y, z))/(kapp_Ez(z, 1)*delta_z);
    
                psi_Hxy(x, y, z) = b_Hy(y, 1)*psi_Hxy(x, y, z) + c_Hy(y, 1)*EzDy;
                psi_Hxz(x, y, z) = b_Hz(z, 1)*psi_Hxz(x, y, z) + c_Hz(z, 1)*EyDz;

                psi_Hyx(x, y, z) = b_Hx(x, 1)*psi_Hyx(x, y, z) + c_Hx(x, 1)*EzDx;
                psi_Hyz(x, y, z) = b_Hz(z, 1)*psi_Hyz(x, y, z) + c_Hz(z, 1)*ExDz;

                psi_Hzx(x, y, z) = b_Hx(x, 1)*psi_Hzx(x, y, z) + c_Hx(x, 1)*EyDx;
                psi_Hzy(x, y, z) = b_Hy(y, 1)*psi_Hzy(x, y, z) + c_Hy(y, 1)*ExDy;
    
                Hx(x, y, z) = Hx(x, y, z) - delta_t/tempmu*( EzDy - EyDz + psi_Hxy(x, y, z) - psi_Hxz(x, y, z) );
                Hy(x, y, z) = Hy(x, y, z) + delta_t/tempmu*( EzDx - ExDz + psi_Hyx(x, y, z) - psi_Hyz(x, y, z) );
                Hz(x, y, z) = Hz(x, y, z) - delta_t/tempmu*( EyDx - ExDy + psi_Hzx(x, y, z) - psi_Hzy(x, y, z) );
            end
        end
    end
   
    if mod(time, plotspeed)==0
        % % change slice of Ez shown
        % % 2D slice
        % imagesc(Ez(:, :, gridsize/2+1));
        % % pcolor(Ez(:, :, gridsize/2+1)');
        % shading flat;
        % view(2);
        % daspect([1,1,1]);
        % colorbar;
        % % clim([-0.4,0.4]);

        slice(Ez,gridsize/2,[],[]);  %slice 2D y-z plot at x-center
        hold on; 
        % slice(Ez,gridsize + 1 - gridsize/4,[],[]);  %slice 2D x-z plot at y-center
        % slice(Ez,[],gridsize/4,[]);  %slice 2D x-y plot at y-center
        % slice(Ez,[],gridsize + 1 - gridsize/4,[]);  %slice 2D x-z plot at y-center
        % slice(Ez,[],[],gridsize/4);  %slice 2D x-y plot at y-center
        % slice(Ez,[],[],gridsize + 1 - gridsize/4);  %slice 2D x-z plot at y-center
        hold off;
        shading flat;
        colorbar;
        daspect([1 1 1]);
        % clim([-0.01 0.01]); 
       

        animation(time/plotspeed + 1, 1) = getframe(f);
    end
    time = time + 1;
end
pause(1);
% animation playback
stop(t);
delete(t);