clc
close all
clear all

%% Initialising 
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

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

%% Calculation
%Exact
[fa, fb] = panelinf(panel_x_min, panel_y_min, ...
    panel_x_max, panel_y_max, xm, ym);
psi_exact = gamma_a .* fa + gamma_b .* fb;

% Approximate with discrete vortices
psi_approx = zeros([ny nx]);

del = norm([panel_x_max - panel_x_min, panel_y_max - panel_y_min]);
nv = 1000;
for i = linspace(0,1,nv)
    gamma = (gamma_a * (1 - i) + gamma_b * i) * del / nv;
    x = panel_x_min * (1 - i) + panel_x_max * i;
    y = panel_y_min * (1 - i) + panel_y_max * i;
    tmp = psipv(x, y, gamma, xm, ym);
    psi_approx = tmp + psi_approx;
end

%% Contour Plots
figure;
contour(xm, ym, psi_exact, c);
title("Contour plot for vortex sheet (exact)", 'FontSize',20);
xlabel("x", 'FontSize',18,'FontWeight','bold');
ylabel("y", 'FontSize',18,'FontWeight','bold');

figure;
contour(xm, ym, psi_approx, c);
title("Contour plot for vortex sheet (approximate)", 'FontSize',20);
xlabel("x", 'FontSize',18,'FontWeight','bold');
ylabel("y", 'FontSize',18,'FontWeight','bold');