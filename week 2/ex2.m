clc
close all
clear all

%% Initial 
Re = 1e6;
n_trapeziums = 100;
x = linspace(0,1,n_trapeziums + 1);
U = 1;
duedx = +0.2;
ue = ones([n_trapeziums + 1,1]);
for i=1:n_trapeziums + 1
    ue(i) = 1 + duedx * (i-1) / n_trapeziums;
end

result = 0;
theta = zeros([n_trapeziums + 1,1]);
laminar = true;
i = 1;
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
        disp([x(i) Rethet/1000]);
    end
end