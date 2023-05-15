clc
close all
clear all

%% Initial 
xmin = -2.5;
xmax = 2.5;
nx = 151;
ymin = -2;
ymax = 2;
ny = 140;
np = 100;

alpha = pi/20;

c = -1.73:0.05:1.75;

R = 1;

%% more initialising

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

theta = (0:np)*2*pi/np;
x_cy = R * cos(theta);
y_cy = R * sin(theta);

%% Streamfunction

A = build_lhs(x_cy,y_cy);
b = build_rhs(x_cy,y_cy,alpha);

gam = A\b;
disp(gam);

psi = ym .* cos(alpha) - xm .* sin(alpha);
for i = 1:(length(x_cy)-1)
    panel_x_min = x_cy(i);
    panel_y_min = y_cy(i);

    panel_x_max = x_cy(i+1);
    panel_y_max = y_cy(i+1);

    [fa, fb] = panelinf(panel_x_min, panel_y_min, ...
    panel_x_max, panel_y_max, xm, ym);
    psi_temp = gam(i) .* fa + gam(i+1) .* fb;
    psi = psi_temp + psi;
end


%% Plotting

figure;
theta_norm = (theta)/pi;

scatter(theta_norm(1:end), gam');

circ = sum(gam)*2*pi*R/np;
disp(circ);

figure;
contour(xm, ym, psi, c);
hold on
plot(x_cy,y_cy)
hold off