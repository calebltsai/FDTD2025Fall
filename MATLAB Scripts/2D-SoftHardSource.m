clear();
w = 2*pi*1e9; % source pulse frequency
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

imagesc(ez_grid);
pcolor(ez_grid');
shading flat;
% surf(ez_grid, 'edgecolor', 'none');
view(2);
daspect([1,1,1]);
colorbar;
clim([-0.2,0.2]);
getframe;

for x=2:400
    for y=2:400
        if ((x-200)**2 + (y-200)**2 <= 20))
            eps

for t=1:600
    %update e with previous e data and 1/2 frame data of h
    
    % soft source
    for x=2:400
        for y=2:400
            ez_grid(x,y) = delta_t/(delta_x*eps)*(hy_grid(x,y)-hy_grid(x-1,y)-hx_grid(x,y)+hx_grid(x,y-1))+ez_grid(x,y);
            if x==200 && t<=400
                ez_grid(x,y) = exp(-((t-100)/50)^2)*exp(-((y-200)/20)^2)*sin(w*(t-100)*delta_t) + ez_grid(x,y);
            end
        end
    end

    % hard source
    % for x=2:400
    %     for y=2:400
    %         if x==200 && t<=400
    %             ez_grid(x,y) = exp(-((t-100)/50)^2)*exp(-((y-200)/20)^2)*sin(w*(t-100)*delta_t);
    %         else
    %             ez_grid(x,y) = delta_t/(delta_x*eps)*(hy_grid(x,y)-hy_grid(x-1,y)-hx_grid(x,y)+hx_grid(x,y-1))+ez_grid(x,y);      
    %         end
    %     end
    % end

    % update h with previous h data and 1/2 frame data of e
    for x=1:400
        for y=1:400
            hy_grid(x,y) = delta_t/(delta_x*mu)*(ez_grid(x+1,y)-ez_grid(x,y))+hy_grid(x,y);

            hx_grid(x,y) = delta_t/(delta_x*mu)*(ez_grid(x,y)-ez_grid(x,y+1))+hx_grid(x,y);

        end
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