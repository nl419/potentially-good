clc
close all
clear all

%% data

matrix1 = xlsread('comparison_camber.xlsx');

index = matrix1(1:7, 1);
cam_pos_var = matrix1(1:7, 12);
LD_thick = matrix1(1:7, 13);
alpha_thick = matrix1(1:7, 14);
cam_var = matrix1(1:7, 2);
LD_cam = matrix1(1:7, 3);
alpha_cam = matrix1(1:7, 4);

disp(index);

%% Plotting

figure
plot(index, LD_cam, '*b', 'linewidth', 1.5);
xlim([0.5 7.5]);
xticklabels(string(cam_var))

xlabel("NACA Aerofoil", 'FontSize',18,'FontWeight','bold');
ylabel("L/D", 'FontSize',18,'FontWeight','bold');

figure
plot(index, alpha_cam, '*r', 'linewidth', 1.5);
xlim([0.5 7.5]);
ylim([0 3]);
xticklabels(string(cam_var))

xlabel("NACA Aerofoil", 'FontSize',18,'FontWeight','bold');
ylabel("$\alpha_{max}$ (degrees)", 'FontSize',18,'FontWeight','bold','interpreter','latex');