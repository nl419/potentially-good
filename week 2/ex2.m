clc
close all
clear all

%% Initial 
Re = 1e5;
n_trapeziums = 5000;
xs_norm = linspace(0,1,n_trapeziums + 1);
U = 1;
gradU = 0.0; % d(ue/U)/d(x/L)
us_norm = ones([n_trapeziums + 1,1]);
for i=1:n_trapeziums + 1
    us_norm(i) = 1 + gradU * (i-1) / n_trapeziums;
end

result = 0;

theta_over_L = zeros([n_trapeziums + 1,1]);

laminar = true;
i = 0;
while laminar && i < n_trapeziums
    i = i + 1;
    result = result + ueintbit(xs_norm(i), us_norm(i), xs_norm(i+1), us_norm(i+1));
    theta_over_L(i+1) = sqrt((0.45 / Re) * power(us_norm(i+1), -6) * result);
    m = -Re * theta_over_L(i+1) * theta_over_L(i+1) * gradU;
    H = thwaites_lookup(m);
    He = laminar_He(H);
    Rethet = Re * us_norm(i+1) * theta_over_L(i+1);

    if log(Rethet) >= 18.4*He - 21.74
        laminar = false;
        disp([xs_norm(i+1) Rethet/1000]);
    end
end

balsius = (0.664 / sqrt(Re)) .* sqrt(xs_norm);


%% Plotting

% disp(balsius);
% disp(theta_over_L);
figure;
scatter(xs_norm, theta_over_L);
hold on
scatter(xs_norm, balsius);
hold off
legend("approximate", "balsius")