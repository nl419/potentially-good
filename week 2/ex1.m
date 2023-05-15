clc
close all
clear all

%% Initial 
Re = 1e5;
n_trapeziums = 100;
xs_norm = linspace(0,1,n_trapeziums + 1);
U = 1;
us_norm = ones([n_trapeziums + 1,1]);

result = 0;

theta_over_L = zeros([n_trapeziums + 1,1]);

for i=1:n_trapeziums
    result = result + ueintbit(xs_norm(i), us_norm(i), xs_norm(i+1), us_norm(i+1));
    theta_over_L(i) = sqrt((0.45 / Re) * power(us_norm(i+1), -6) * result);
end

balsius = (0.664 / sqrt(Re)) .* sqrt(xs_norm);


%% Plotting

disp(balsius);
disp(theta_over_L);
figure;
scatter(xs_norm, theta_over_L);
hold on
scatter(xs_norm, balsius);
hold off