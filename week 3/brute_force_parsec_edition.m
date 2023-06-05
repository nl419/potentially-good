clc
close all
clear all

Re = 0.5e6;
np = 400;
alpha = 0:0.5:10;

rle = linspace(0.001, 0.03, 10);
x_low = linspace(0, 1, 10);
y_low = linspace(-0.2, 0.2, 10);
d2ydx2_low = linspace(-0.4, 0.4, 10);
dth_low = param(5);
dx_up = param(6);
dy_up = param(7);
dd2ydx2_up = param(8);
dth_up = param(9);

cambers = 0:9;
thicks = 1:30;
camber_poss = 1:9;

ids = strings([length(cambers), length(camber_poss), length(thicks)]);

for icamber = 1:length(cambers)
    for icamber_pos = 1:length(camber_poss)
        for ithick = 1:length(thicks)
            ids(icamber, icamber_pos, ithick) = strcat(num2str(cambers(icamber), "%01.f"), ...
                num2str(camber_poss(icamber_pos), "%01.f"), ...
                num2str(thicks(ithick), "%02.f"));
        end
    end
end
ids = reshape(ids, [], 1);

L_over_Ds = zeros(size(ids));

tic
parfor i = 1:length(L_over_Ds)
    L_over_D = foil_L_over_D(char(ids(i)), np, alpha, Re);
    L_over_Ds(i) = L_over_D;
end
clc
toc
best_L_over_D = max(L_over_Ds);
best_ids = ids(L_over_Ds==best_L_over_D);
disp(best_L_over_D);
disp(best_ids);

