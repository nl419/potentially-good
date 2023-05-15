clc
close all
clear all

%% Initialising
xmin = -2.5;
xmax = 2.5;
nx = 51;
ymin = -2;
ymax = 2;
ny = 42;
gamma_a = 3.0;
gamma_b = 4.0;

del = 1.5;
c = -1:0.05:1;

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

%% Calculation
%Exact solution
[fa, fb] = refpaninf(del, xm, ym);
psi_exact = gamma_a .* fa + gamma_b .* fb;

%Approximate with discrete vortices
psi_approx = zeros([ny nx]);
panel_x_min = 0;
panel_x_max = del + panel_x_min;
nv = 100;
for i = linspace(0,1,nv)
    gamma = (gamma_a * (1 - i) + gamma_b * i) * del / nv;
    x = panel_x_min * (1 - i) + panel_x_max * i;
    tmp = psipv(x, 0, gamma, xm, ym);
    psi_approx = tmp + psi_approx;
end

%% Plotting
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