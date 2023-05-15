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

xc = 0.75;
yc = 0.5;

Gamma = 3.0;

c = -0.4:0.2:1.2;

%% Actually doing stuff

x_1D = linspace(xmin, xmax, nx);
y_1D = linspace(ymin, ymax, ny);
xm = ones([ny 1])*x_1D;
ym = y_1D'*ones([1 nx]);

psi = psipv(xc, yc, Gamma, xm, ym);

% asdf = reshape(xm, [], 1);
% disp(asdf)

contour(xm,ym,psi,c);