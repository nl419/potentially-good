clc
close all
clear all

%% Initial

global Re ue0 duedx

Re = 1e6;
ue0 = 1;
duedx = 0;

x0 = 0.01;
thick0(1) = 0.037 * x0 * (Re * x0) ^ (-1/5);
thick0(2) = 1.80 * thick0(1);

[delx, thickhist] = ode45(@thickdash, [0 0.99], thick0);
disp(delx);
disp(thickhist);

x = x0 + delx;
theta_over_L = thickhist(:,1);

theta_7 = 0.037*x.*(Re.*x).^(-1/5);
theta_9 = 0.023*x.*(Re.*x).^(-1/6);

%% Plotting

figure;
scatter(x, theta_over_L);
hold on
scatter(x, theta_7);
hold on
scatter(x, theta_9);

legend('Approx', 'theta 7', 'theta 9');
