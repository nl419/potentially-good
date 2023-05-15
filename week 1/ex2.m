clc
close all
clear all

%% Initial 
xmin = -2.5;
xmax = 2.5;
nx = 51;
ymin = -2;
ymax = 2;
ny = 41;

gamma_a = 3.0;
gamma_b = 4.0;

del = 1.5;

c = -1:0.05:1; % TODO check this 

%% more initialising

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

%% Exact solution
[fa, fb] = refpaninf(del, xm, ym);
psi_good = gamma_a .* fa + gamma_b .* fb;

%% Approximate with discrete vortices
psi_bad = zeros([ny nx]);
panel_x_min = 0;
panel_x_max = del + panel_x_min;
nv = 100;
for i = linspace(0,1,nv)
    gamma = (gamma_a * (1 - i) + gamma_b * i) * del / nv;
    x = panel_x_min * (1 - i) + panel_x_max * i;
    tmp = psipv(x, 0, gamma, xm, ym);
    psi_bad = tmp + psi_bad;
end

%% looking at the horrors we have produced
contour(xm, ym, psi_good, c);
figure;
contour(xm, ym, psi_bad, c);