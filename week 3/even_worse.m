clc
close all
clear all

%% Init

Re = 0.5e6;
np = 400;
alpha = 0:0.25:8;
id = '0020';

%% Literal hell
[xk0, yk0] = naca4_nowrite(id);
% [xk0, yk0] = textread("Geometry\opt.surf", '%f%f');

uplim_y = 0.05*ones([length(yk0),1]);
lowlim_y = -0.05*ones([length(yk0),1]);
uplim_y(1) = 0;
uplim_y(end) = 0;
lowlim_y(1) = 0;
lowlim_y(end) = 0;
uplim_y = uplim_y + yk0;
lowlim_y = lowlim_y + yk0;

opts = optimoptions(@fmincon,'MaxFunctionEvaluations',500,'UseParallel',true);

y_fin = fmincon(@(y)foil_L_over_D_noid_negative(xk0, y, np, alpha, Re, true),yk0, ...
    [],[],[],[],lowlim_y,uplim_y,[],opts);
max_LD = -foil_L_over_D_noid_negative(xk0, y_fin, np, alpha, Re, false);

disp(max_LD);
disp(-foil_L_over_D_noid_negative(xk0, yk0, np, alpha, Re, false));
plot(xk0, y_fin);
hold on;
plot(xk0, yk0);
legend("New", "Old");

fname = ['Geometry/opt_brick.surf'];

fid = fopen(fname,'w');
fprintf(fid, '%10.6f %10.6f \n', [xk0';y_fin']);
fclose(fid);
