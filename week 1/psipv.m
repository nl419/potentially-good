function [psi] = psipv(xc, yc, Gamma, xs, ys)
x2 = xs - xc;
y2 = ys - yc;
r_squared = x2.*x2 + y2.*y2;
psi = (-Gamma / (4 * pi)) * log(r_squared);
end