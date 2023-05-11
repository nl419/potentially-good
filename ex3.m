clc
close all
clear all

%% Initial 
xmin = 0;
xmax = 5;
nx = 151;
ymin = 0;
ymax = 4;
ny = 140;

gamma_a = 3.0;
gamma_b = 4.0;

del = 1.5;

c = -1:0.05:1;

panel_x_min = 2.2;
panel_x_max = 4.1;
panel_y_min = 1.3;
panel_y_max = 2.9;

%% more initialising

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

%% Exact solution
[fa, fb] = panelinf(panel_x_min, panel_y_min, ...
    panel_x_max, panel_y_max, xm, ym);
psi_good = gamma_a .* fa + gamma_b .* fb;

%% Approximate with discrete vortices
psi_bad = zeros([ny nx]);

del = norm([panel_x_max - panel_x_min, panel_y_max - panel_y_min]);
nv = 1000;
for i = linspace(0,1,nv)
    gamma = (gamma_a * (1 - i) + gamma_b * i) * del / nv;
    x = panel_x_min * (1 - i) + panel_x_max * i;
    y = panel_y_min * (1 - i) + panel_y_max * i;
    tmp = psipv(x, y, gamma, xm, ym);
    psi_bad = tmp + psi_bad;
end

%% looking at the horrors we have produced
contour(xm, ym, psi_good, c);
figure;
contour(xm, ym, psi_bad, c);