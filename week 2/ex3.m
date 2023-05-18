clc
close all
clear all

%% Initial 
Re = 1.83e6;
n_trapeziums = 5000;
xs_norm = linspace(0,1,n_trapeziums + 1);
gradU_norm = -0.5; % d(ue/U)/d(x/L)
us_norm = ones([n_trapeziums + 1,1]);
for i=1:n_trapeziums + 1
    us_norm(i) = 1 + gradU_norm * (i-1) / n_trapeziums;
end

result = 0;

theta_over_L = zeros([n_trapeziums + 1,1]);

laminar = true;
i = 0;

int = 0; 
ils = 0;
Rethet = 0;

while laminar && i < n_trapeziums
    i = i + 1;
    result = result + ueintbit(xs_norm(i), us_norm(i), xs_norm(i+1), us_norm(i+1));
    theta_over_L(i+1) = sqrt((0.45 / Re) * power(us_norm(i+1), -6) * result);
    m = -Re * theta_over_L(i+1) * theta_over_L(i+1) * gradU_norm;
    H = thwaites_lookup(m);
    He = laminar_He(H);
    Rethet = Re * us_norm(i+1) * theta_over_L(i+1);

    if log(Rethet) >= 18.4*He - 21.74
        laminar = false;
        int = xs_norm(i+1);
    elseif m >= 0.09
        laminar = false;
        ils = xs_norm(i+1);
    end
end

balsius = (0.664 / sqrt(Re)) .* sqrt(xs_norm);

if int ~= 0
    disp(['Natural transition at ' num2str(int) ' with Rethet ' num2str(Rethet)]);
end
if ils ~= 0
    disp(['Laminar separation at ' num2str(ils) ' with Rethet ' num2str(Rethet)]);
end

%% Plotting

% disp(balsius);
% disp(theta_over_L);
% figure;
% scatter(xs_norm, theta_over_L);
% hold on
% scatter(xs_norm, balsius);
% hold off
% legend("approximate", "balsius")