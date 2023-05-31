clc
close all
clear all

%% Init

Re = 0.5e6;
np = 400;
alpha = 0:1:20;
id = '8505';

%% Literal hell
[xk0, yk0] = naca4_nowrite(id);

uplim_y = 0.15*ones(length(yk0));
lowlim_y = -0.05*ones(length(yk0));

fun = @foil_L_over_D_noid;
y_fin = fmincon(@(y)foil_L_over_D_noid(xk0, y, np, alpha, Re),yk0,[],[],[],[],lowlim_y,uplim_y);
max_LD = foil_L_over_D_noid(xk0, y_fin, np, alpha, Re);

disp(-max_LD);