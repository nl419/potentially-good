function f = ueintbit(xa_norm,ua_norm,xb_norm,ub_norm)
    ubar = (ua_norm + ub_norm) ./ 2;
    delu = ub_norm - ua_norm;
    delx = xb_norm - xa_norm;
    f = power(ubar, 5) + (5 / 6) .* power(ubar, 3) .* power(delu, 2) ...
            + (1 / 16) * ubar * power(delu, 4);
    f = f * delx;
end