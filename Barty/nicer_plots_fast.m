clc
close all
clear all

%% data

matrix1 = xlsread('comparison_highspeed.xlsx');

index = matrix1(1:7, 1);

cam_var = matrix1(1:7, 2);
LD_cam = matrix1(1:7, 3);
alpha_cam = matrix1(1:7, 4);

cam_pos_var = matrix1(1:7, 7);
LD_cam_pos = matrix1(1:7, 8);
alpha_cam_pos = matrix1(1:7, 9);

thick_var = matrix1(1:7, 12);
LD_thick = matrix1(1:7, 13);
alpha_thick = matrix1(1:7, 14);


disp(index);

%% Plotting

figure
plot(index, LD_thick, '*b', 'linewidth', 1.5);
xlim([0.5 7.5]);
xticklabels(string(thick_var))

xlabel("NACA Aerofoil", 'FontSize',18,'FontWeight','bold');
ylabel("L/D", 'FontSize',18,'FontWeight','bold');

figure
plot(index, alpha_thick, '*r', 'linewidth', 1.5);
xlim([0.5 7.5]);
ylim([0 9]);
xticklabels(string(thick_var))

xlabel("NACA Aerofoil", 'FontSize',18,'FontWeight','bold');
ylabel("$\alpha_{max}$ (degrees)", 'FontSize',18,'FontWeight','bold','interpreter','latex');