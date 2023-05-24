clc
close all
clear all

%% Initial 
Re = 1e5;
n_trapeziums = 100;
x = linspace(0,1,n_trapeziums + 1);
duedx = -0.5; % d(ue/U)/d(x/L)
ue = ones([n_trapeziums + 1,1]);
for i=1:n_trapeziums + 1
    ue(i) = 1 + duedx * (i-1) / n_trapeziums;
end

result = 0;

theta = zeros([n_trapeziums + 1,1]);

laminar = true;
i = 1;

int = 0; 
ils = 0;
Rethet = 0;

while laminar && i < (n_trapeziums+1)
    i = i + 1;
    result = result + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta(i) = sqrt((0.45 / Re) * power(ue(i), -6) * result);
    m = -Re * theta(i) * theta(i) * duedx;
    H = thwaites_lookup(m);
    He = laminar_He(H);
    Rethet = Re * ue(i) * theta(i);

    if log(Rethet) >= 18.4*He - 21.74
        laminar = false;
        int = x(i);
    elseif m >= 0.09
        laminar = false;
        ils = x(i);
    end
end

if int ~= 0
    disp(['Natural transition at ' num2str(int) ' with Rethet ' num2str(Rethet)]);
end
if ils ~= 0
    disp(['Laminar separation at ' num2str(ils) ' with Rethet ' num2str(Rethet)]);
end