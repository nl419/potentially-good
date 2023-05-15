clc
close all
clear all

%% Initialisation 
xmin = -2.5;
xmax = 2.5;
nx = 51;
ymin = -2;
ymax = 2;
ny = 41;
xc = 0.75;
yc = 0.5;

Gamma = 3.0;
c = -0.4:0.2:1.2;

%% Calculation

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);
psi = psipv(xc, yc, Gamma, xm, ym);

%% Plotting

contour(xm,ym,psi,c);
title("Contour plot of a point vortex", 'FontSize',20);
xlabel("x", 'FontSize',18,'FontWeight','bold');
ylabel("y", 'FontSize',18,'FontWeight','bold');