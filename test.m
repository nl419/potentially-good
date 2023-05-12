clc
close all
clear all

np =100;
theta = (0:np)*2*pi/np;
x_cy = cos(theta);
y_cy = sin(theta);
alpha = 0;

A = [1 2 3 4;
    5 6 7 8;
    9 10 11 12];

disp(size(A));

psi_n2 = circshift(A,-1,1);
disp(psi_n2);


