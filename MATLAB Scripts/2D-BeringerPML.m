clear();
w = 2*pi*1e9;
c = 299792458;
lambda = 2*pi*c/w;
mu = 4*pi*1e-7;
eps = 1/(c^2*mu);
delta_x = lambda/50;
delta_t = 0.5*delta_x/c;

gridsize = 100;
plotspeed = 3;
tolerance = 1e-15;
m = 2.5;
pmlwidth = 10;

%initialize with blank frame
ez_grid = zeros(gridsize + 1);
ez_x = zeros(gridsize + 1);
ez_y = zeros(gridsize + 1);

hy_grid = zeros(gridsize);
hx_grid = zeros(gridsize);

% testing reflective error of pml
large_ez = zeros(gridsize*3+1);
large_hy = zeros(gridsize*3);
large_hx = zeros(gridsize*3);

% build array of coeffs
eta = sqrt(mu/eps);
sig_e = -log(tolerance)*(m+1)/(2*eta*pmlwidth*delta_x);
sig_m = sig_e*eta*eta;

ca_x = ones(gridsize+1, 1)*delta_t/eps;
cb_x = ones(gridsize+1, 1)*eps/delta_t;
ca_y = ones(gridsize+1, 1)*delta_t/eps;
cb_y = ones(gridsize+1, 1)*eps/delta_t;

da_x = ones(gridsize, 1)*delta_t/mu;
db_x = ones(gridsize, 1)*mu/delta_t;
da_y = ones(gridsize, 1)*delta_t/mu;
db_y = ones(gridsize, 1)*mu/delta_t;


% arrays for pml boundary resable in square grid and uniform boundary
% widths, coeffs increase corresponding to index
% (x)/pmlwidth
for x=1:pmlwidth
    ca_x(pmlwidth + 2 - x, 1) = 1/(eps/delta_t+(sig_e*((x - 1)/pmlwidth).^m)/2);
    ca_x(gridsize - pmlwidth + x, 1) = 1/(eps/delta_t+(sig_e*((x - 1)/pmlwidth).^m)/2); % change index
    cb_x(pmlwidth + 2 - x, 1) = (eps/delta_t-(sig_e*((x - 1)/pmlwidth).^m)/2);
    cb_x(gridsize - pmlwidth + x, 1) = (eps/delta_t-(sig_e*((x - 1)/pmlwidth).^m)/2); % change index

    ca_y(pmlwidth + 2 - x, 1) = 1/(eps/delta_t+(sig_e*((x - 1)/pmlwidth).^m)/2);
    ca_y(gridsize - pmlwidth + x, 1) = 1/(eps/delta_t+(sig_e*((x - 1)/pmlwidth).^m)/2); % change index 
    cb_y(pmlwidth + 2 - x, 1) = (eps/delta_t-(sig_e*((x - 1)/pmlwidth).^m)/2);
    cb_y(gridsize - pmlwidth + x, 1) = (eps/delta_t-(sig_e*((x - 1)/pmlwidth).^m)/2); % change index
        
    da_x(pmlwidth + 1 - x, 1) = 1/(mu/delta_t+sig_m*((x - 0.5)/pmlwidth).^m/2);
    da_x(gridsize - pmlwidth + x, 1) = 1/(mu/delta_t+sig_m*((x - 0.5)/pmlwidth).^m/2);
    da_y(pmlwidth + 1 - x, 1) = 1/(mu/delta_t+sig_m*((x - 0.5)/pmlwidth).^m/2); % change index
    da_y(gridsize - pmlwidth + x, 1) = 1/(mu/delta_t+sig_m*((x - 0.5)/pmlwidth).^m/2);
    
    db_x(pmlwidth + 1 - x, 1) = (mu/delta_t-sig_m*((x - 0.5)/pmlwidth).^m/2); 
    db_x(gridsize - pmlwidth + x, 1) = (mu/delta_t-sig_m*((x - 0.5)/pmlwidth).^m/2); 
    db_y(pmlwidth + 1 - x, 1) = (mu/delta_t-sig_m*((x - 0.5)/pmlwidth).^m/2); 
    db_y(gridsize - pmlwidth + x, 1) = (mu/delta_t-sig_m*((x - 0.5)/pmlwidth).^m/2); 
end

% initial frame
imagesc(ez_grid);
pcolor(ez_grid');
shading flat;
% surf(ez_grid, 'edgecolor', 'none');
view(2);
daspect([1,1,1]);
colorbar;
clim([-0.2,0.2]);
getframe;
pause(0.1*plotspeed)


for t=1:300
    %update e with previous e data and 1/2 frame data of h
    
    % soft source
    for x=2:gridsize
        for y=2:gridsize
            ez_x(x,y) = ca_y(y, 1)*(cb_y(y, 1)*ez_x(x,y)+(hx_grid(x,y-1)-hx_grid(x,y))/delta_x);
            ez_y(x,y) = ca_x(x, 1)*(cb_x(x, 1)*ez_y(x,y)+(hy_grid(x,y)-hy_grid(x-1,y))/delta_x);
            ez_grid(x,y) =  ez_y(x,y) + ez_x(x,y);
        end
    end
    ez_grid(gridsize/2,gridsize/2) = exp(-((t-100)/30)^2)*sin(w*(t-100)*delta_t) + ez_grid(gridsize,gridsize);

    % update h with previous h data and 1/2 frame data of e
    for x=1:gridsize
        for y=1:gridsize
            hy_grid(x,y) = da_x(x, 1)*(db_x(x, 1)*hy_grid(x,y)+(ez_grid(x+1,y)-ez_grid(x,y))/delta_x);
            hx_grid(x,y) = da_y(y, 1)*(db_y(y, 1)*hx_grid(x,y)+(ez_grid(x,y)-ez_grid(x,y+1))/delta_x);
        end
    end


    for x=2:gridsize*3
        for y=2:gridsize*3
            large_ez(x,y) = delta_t/(delta_x*eps)*(large_hy(x,y)-large_hy(x-1,y)-large_hx(x,y)+large_hx(x,y-1))+large_ez(x,y);
        end
    end
    
    large_ez(gridsize*3/2,gridsize*3/2) = large_ez(gridsize*3/2, gridsize*3/2) + exp(-((t-100)/30)^2)*sin(w*(t-100)*delta_t);
    
    for x=1:gridsize*3
        for y=1:gridsize*3
            large_hy(x,y) = large_hy(x,y)+delta_t/(delta_x*mu)*(large_ez(x+1,y)-large_ez(x,y));
            large_hx(x,y) = large_hx(x,y)+delta_t/(delta_x*mu)*(large_ez(x,y)-large_ez(x,y+1));
        end
    end
    
    submatrix = large_ez(gridsize+1:gridsize*2+1, gridsize+1:gridsize*2+1);
   
    if mod(t, plotspeed)==0
        imagesc(ez_grid);
        pcolor(ez_grid');
        % imagesc(large_ez);
        % pcolor(large_ez');
        shading flat;
        % surf(ez_grid, 'edgecolor', 'none');
        view(2);
        daspect([1,1,1]);
        colorbar;
        clim([-0.2,0.2]);
        
        getframe;
        pause(0.1*plotspeed)
    end
end