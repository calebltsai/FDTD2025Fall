clear();
w = 2*pi*1e9;
c = 299792458;
lambda = 2*pi*c/w;
mu = 4*pi*1e-7;
eps = 1/(c^2*mu);
delta_x = lambda/40;
delta_t = delta_x/c;

%time_grid = zeros(4000,4000,2);
%or graph each e frame and save for movie

% scattered field at egrid 301 and 100 respective h scattered to left at 100 and
% right at 302 and 

%initialize with blank frame
ez_grid = zeros(401,1);
hy_grid = zeros(400,1);

inc_ez = zeros(401,1);
inc_hy = zeros(400,1);

plot(ez_grid);
xlim([0 402]);
ylim([-1 1]);
animation = [getframe];
for t=1:1000
    %update e with previous e data and 1/2 frame data of h
    % ez_grid(1,1) = sin(w*t*delta_t);
    
    for x=2:400
        ez_grid(x,1) = delta_t/(delta_x*eps)*(hy_grid(x,1)-hy_grid(x-1,1))+ez_grid(x,1);
        inc_ez(x,1) = delta_t/(delta_x*eps)*(inc_hy(x,1)-inc_hy(x-1,1))+inc_ez(x,1);
    end
    ez_grid(200,1) = ez_grid(200,1) + exp(-(50-t)^2/90)*sin(w*(50-t)*delta_t);
    inc_ez(200,1) = inc_ez(200,1) + exp(-(50-t)^2/90)*sin(w*(50-t)*delta_t);

    ez_grid(250,1) = 0;

    ez_grid(301,1) = ez_grid(301, 1) + delta_t/(delta_x*eps)*inc_hy(300,1);
    ez_grid(100,1) = ez_grid(100, 1) - delta_t/(delta_x*eps)*inc_hy(100,1);
        
    %update h with previous h data and 1/2 frame data of e
    for x=1:400
        hy_grid(x,1) = delta_t/(delta_x*mu)*(ez_grid(x+1,1)-ez_grid(x,1))+hy_grid(x,1);
        inc_hy(x,1) = delta_t/(delta_x*mu)*(inc_ez(x+1,1)-inc_ez(x,1))+inc_hy(x,1);
    end
    hy_grid(300,1) = hy_grid(300, 1) + delta_t/(delta_x*mu)*inc_ez(301,1);
    hy_grid(100,1) = hy_grid(100, 1) - delta_t/(delta_x*mu)*inc_ez(100,1);

    plot(ez_grid);
    xlim([0 402]);
    ylim([-1 1]);
    animation(t,1) = getframe;
end
movie(animation, 1, 40);