clc
close all
clear all

%% Initialising

Re = 1e5;
n_trapeziums = 100;
xs_norm = linspace(0,1,n_trapeziums + 1);
U = 1;
us_norm = ones([n_trapeziums + 1,1]);
result = 0;
theta_over_L = zeros([n_trapeziums + 1,1]);
theta_over_L(1) = result;
for i=1:n_trapeziums
    result = result + ueintbit(xs_norm(i), us_norm(i), xs_norm(i+1), us_norm(i+1));
    theta_over_L(i+1) = sqrt((0.45 / Re) * power(us_norm(i+1), -6) * result);
end
blasius = (0.664 / sqrt(Re)) .* sqrt(xs_norm);

%% Plotting

figure;
scatter(xs_norm, theta_over_L);
hold on
scatter(xs_norm, blasius);
hold off
legend("Approximate Solution", "Blasius Solution")
xlabel("x/L", 'FontSize',18,'FontWeight','bold');
ylabel("$\theta$", 'interpreter','latex', 'FontSize',18,'FontWeight','bold');