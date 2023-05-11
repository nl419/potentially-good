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

c = -1.75:0.05:1.75;

%% more initialising

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

theta = (0:np)*2*pi/np;
x_cy = cos(theta);
y_cy = sin(theta);

gamma = -2*sin(theta);

%% Streamfunction

psi = ym;
for i = 1:(length(x_cy)-1)
    panel_x_min = x_cy(i);
    panel_y_min = y_cy(i);

    panel_x_max = x_cy(i+1);
    panel_y_max = y_cy(i+1);

    [fa, fb] = panelinf(panel_x_min, panel_y_min, ...
    panel_x_max, panel_y_max, xm, ym);
    psi_temp = gamma(i) .* fa + gamma(i+1) .* fb;
    psi = psi_temp + psi;
end

%% Plotting
figure;
contour(xm, ym, psi, c);
hold on
plot(x_cy,y_cy)
hold off