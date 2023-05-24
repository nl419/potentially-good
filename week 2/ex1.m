clc
close all
clear all

%% Initialising

Re = 2500;
n_trapeziums = 100;
x = linspace(0,1,n_trapeziums + 1);
U = 1;
ue = ones([n_trapeziums + 1,1]);
result = 0;
theta = zeros([n_trapeziums + 1,1]);
theta(1) = result;
for i=1:n_trapeziums
    result = result + ueintbit(x(i),ue(i), x(i+1),ue(i+1));
    theta(i+1) = sqrt((0.45 / Re) * power(ue(i+1), -6) * result);
end
blasius = (0.664 / sqrt(Re)) .* sqrt(x);

%% Plotting

figure;
scatter(x, theta);
hold on
scatter(x, blasius);
hold off
legend("Approximate Solution", "Blasius Solution",'FontSize',14)
xlabel("x/L", 'FontSize',18,'FontWeight','bold');
ylabel("$\theta/L$", 'interpreter','latex', 'FontSize',18,'FontWeight','bold');