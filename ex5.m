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

alpha = 0.1;

c = -1.75:0.05:1.75;

%% more initialising

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

theta = (0:np)*2*pi/np;
x_cy = cos(theta);
y_cy = sin(theta);

%% Streamfunction

A = build_lhs(x_cy,y_cy);
b = build_rhs(x_cy,y_cy,alpha);

gam = A\b;
disp(gam);

%% Plotting
