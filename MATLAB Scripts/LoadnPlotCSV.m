clear;

TE = readmatrix("FDTD method\Data\2DTESpectrum.csv");

figure(1);
surf(TE', 'edgecolor', 'none');
% pcolor(TE');
% shading flat;
ylim([1 40]);
% 
% clim([0 2*mean(max(TE(:,1:30)))]);

% TM = readmatrix("FDTD method\Data\2DTMSpectrum.csv");

% figure(2);
% pcolor(TM');
% shading flat;
% ylim([0 500]);
% 
% clim([0 max(mean(TM(:,1:500)))]);

D = readmatrix("FDTD method\Data\dielectrics.csv");

figure(2);
pcolor(D');
shading flat;
axis equal;
