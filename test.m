clc
close all
clear all

np =100;
theta = (0:np)*2*pi/np;
x_cy = cos(theta);
y_cy = sin(theta);
alpha = 0;


rhsvec = build_rhs(x_cy,y_cy,alpha);
disp(rhsvec);