clear();
w = 2*pi*1e9;
c = 299792458;
lambda = 2*pi*c/w;
mu = 4*pi*1e-7;
eps = 1/(c^2*mu);
delta_x = lambda/50;
delta_t = 0.5*delta_x/c;


%initialize with blank frame
ez_grid = zeros(401);
hy_grid = zeros(400);
hx_grid = zeros(400);

inc_ez_grid = zeros(401,1);
inc2 = 0;
inc2_old = 0;
inc2right = 0;
inc2_oldright = 0;
inc_hy_grid = zeros(400,1);


imagesc(ez_grid);
pcolor(ez_grid');
shading flat;
% surf(ez_grid, 'edgecolor', 'none');
view(2);
daspect([1,1,1]);
colorbar;
clim([-0.2,0.2]);
getframe;
for t=1:800
    %update e with previous e data and 1/2 frame data of h
    % soft source
    for x=2:400
        inc_ez_grid(x,1) = inc_ez_grid(x,1) + delta_t/(delta_x*eps)*(inc_hy_grid(x,1)-inc_hy_grid(x-1,1));
        for y=2:400
            % update incident wave grid independent of 
            
            % update visualized grid, soft source
            % box pec
            ez_grid(x,y) = delta_t/(delta_x*eps)*(hy_grid(x,y)-hy_grid(x-1,y)-hx_grid(x,y)+hx_grid(x,y-1))+ez_grid(x,y);
        end
    end
    inc_ez_grid(1,1) = inc2_old;
    inc2_old = inc2;
    inc2 = inc_ez_grid(2,1);
    
    inc_ez_grid(401,1) = inc2_oldright;
    inc2_oldright = inc2right;
    inc2right = inc_ez_grid(400,1);

    for j=125:275
        ez_grid(275,j) = ez_grid(275,j) + delta_t/(delta_x*eps)*(inc_hy_grid(275,1));
        ez_grid(125,j) = ez_grid(125,j) - delta_t/(delta_x*eps)*(inc_hy_grid(124,1));
    end
    inc_ez_grid(100,1) = exp(-((t-100)/50)^2)*sin(w*(t-100)*delta_t) + inc_ez_grid(100,1);
    
    ez_grid(200:210,200:210) = 0;

    % update h with previous h data and 1/2 frame data of e
    for x=1:400
        inc_hy_grid(x,1) = delta_t/(delta_x*mu)*(inc_ez_grid(x+1,1)-inc_ez_grid(x,1))+inc_hy_grid(x,1);
        for y=1:400
            hy_grid(x,y) = delta_t/(delta_x*mu)*(ez_grid(x+1,y)-ez_grid(x,y))+hy_grid(x,y);
            hx_grid(x,y) = delta_t/(delta_x*mu)*(ez_grid(x,y)-ez_grid(x,y+1))+hx_grid(x,y);
        end
    end

    % bounday corrections in transition region
    % hy scattered field at x = 125, 275
    % hx scattered field at y = 125, 275
    for j=125:275
        hy_grid(275,j) = hy_grid(275,j) + delta_t/(delta_x*mu)*(inc_ez_grid(275,1));
        hy_grid(124,j) = hy_grid(124,j) - delta_t/(delta_x*mu)*(inc_ez_grid(125,1));
    end
    for i=125:275
        hx_grid(i,275) = hx_grid(i,275) - delta_t/(delta_x*mu)*(inc_ez_grid(i,1));
        hx_grid(i,124) = hx_grid(i,124) + delta_t/(delta_x*mu)*(inc_ez_grid(i,1));
    end

    if mod(t,10)==0   
        imagesc(ez_grid);
        pcolor(ez_grid');
        shading flat;
        % surf(ez_grid, 'edgecolor', 'none');
        view(2);
        daspect([1,1,1]);
        colorbar;
        clim([-0.2,0.2]);
        getframe;
    end
end