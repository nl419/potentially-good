clc
close all
clear all

%% Initial

global Re ue0 duedx

Re = 1e7;
ue0 = 1;
duedx = -0.9;

x0 = 0.01;
thick0(1) = 0.037 * x0 * (Re * x0) ^ (-1/5);
thick0(2) = 1.80 * thick0(1);

[delx, thickhist] = ode45(@thickdash, [0 0.99], thick0);

x = x0 + delx;
theta_over_L = thickhist(:,1);
delta_over_L = thickhist(:,2);

theta_7 = 0.037*x.*(Re.*x).^(-1/5);
theta_9 = 0.023*x.*(Re.*x).^(-1/6);

He = theta_over_L ./ delta_over_L;
for i = 1:length(He)
   if He(i) >= 1.46
      disp([x(i), He(i)])
      break
   end
end


