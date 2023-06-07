clc
close all
clear all

%% data

load("naca7504.mat");
alpha_old = alpha;
cpu_old = cpu_max;
cpl_old = cpl_max;
cd_old = cdswp;
cl_old = clswp;
LD_old = lovdswp;
su_old = su_max;
sl_old = sl_max;
xs_old = xs;
ys_old = ys;

load("barty5.mat");
alpha_new = alpha;
cpu_new = cpu_max;
cpl_new = cpl_max;
cd_new = cdswp;
cl_new = clswp;
LD_new = lovdswp;
su_new = su_max;
sl_new = sl_max;
xs_new = xs;
ys_new = ys;

%% Plotting

figure;
plot(alpha_new(1:24), LD_new(1:24), '-*b', 'linewidth', 1.5);
hold on
plot(alpha_old(1:24), LD_old(1:24), '-*r', 'linewidth', 1.5);
hold off
xlabel("$\alpha$ (degrees)", 'FontSize',18,'FontWeight','bold','interpreter','latex');
ylabel("L/D", 'FontSize',18,'FontWeight','bold');
legend('Modified Aerofoil', 'NACA 7504');

figure;
plot(su_new, -cpu_new, '-b', 'linewidth', 1.5);
hold on
plot(su_old, -cpu_old, '-r', 'linewidth', 1.5);
hold off
xlabel("$x$", 'FontSize',18,'FontWeight','bold','interpreter','latex');
ylim([-0.2 1]);
ylabel("$-cp_u$", 'FontSize',18,'FontWeight','bold','interpreter','latex');
legend('Modified Aerofoil', 'NACA 7504');

figure;
plot(sl_new, -cpl_new, '-b', 'linewidth', 1.5);
hold on
plot(sl_old, -cpl_old, '-r', 'linewidth', 1.5);
hold off
xlabel("$x$", 'FontSize',18,'FontWeight','bold','interpreter','latex');
ylim([-0.5 -0.1]);
ylabel("$-cp_l$", 'FontSize',18,'FontWeight','bold','interpreter','latex');
legend('Modified Aerofoil', 'NACA 7504');

figure
plot(xs_new, ys_new, '-b', 'linewidth', 1.5);
hold on
plot(xs_old, ys_old, '-r', 'linewidth', 1.5);
hold off
xlabel("$x$", 'FontSize',18,'FontWeight','bold','interpreter','latex');
ylabel("$y$", 'FontSize',18,'FontWeight','bold','interpreter','latex');
legend('Modified Aerofoil', 'NACA 7504');
