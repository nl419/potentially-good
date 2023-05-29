function dthickdx = thickdash(xmx0, thick, Re, ue0, duedx)
    He = thick(2) / thick(1);
    if He >= 1.46
        H = (11 * He + 15) / (48 * He - 59);
    else
        H = 2.803;
    end
    ue = ue0 + duedx * xmx0;
    Rethet = Re * ue * thick(1);
    cf = 0.091448 * power((H - 1) * Rethet, -0.232) * exp(-1.26 * H);
    cdiss = 0.010023 * power((H - 1) * Rethet, -1/6);
    f1 = cf / 2 - (H + 2) * duedx * thick(1) / ue;
    f2 = cdiss - 3 * duedx * thick(2) / ue;
    dthickdx = [f1 f2]';
end