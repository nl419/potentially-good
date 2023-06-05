clc
close all
clear all

%% Init

Re = 0.5e6;
np = 400;
alpha = 0:0.25:8;

%% Literal hell
opts = optimoptions(@fmincon,'MaxFunctionEvaluations',500,'UseParallel',false);

param0 = [0.01,0.450,-0.006,-0.2,0.05,0.350,0.055,-0.350,-6];
param_lowlim = [0, 0, -1, -10, -100, 0, -1, -10, -100];
param_uplim = [100, 1, 1,  10,  100, 1,  1,  10,  100];

% param_fin = fmincon(@(param)foil_L_over_D_parsec(param, np, alpha, Re),param0, ...
%     [],[],[],[],[],[],[],opts);
param_fin = fmincon(@(param)foil_L_over_D_parsec(param, np, alpha, Re),param0, ...
    [],[],[],[],param_lowlim,param_uplim,[],opts);
max_LD = -foil_L_over_D_parsec(param_fin, np, alpha, Re);

disp(max_LD);
disp(-foil_L_over_D_parsec(param0, np, alpha, Re));
% plot(xk0, y_fin);
% hold on;
% plot(xk0, yk0);
% legend("New", "Old");
