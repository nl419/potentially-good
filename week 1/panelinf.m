function [infa, infb] = panelinf(xa, ya, xb, yb, xm, ym)
    t = [(xb - xa) (yb - ya)];
    del = norm(t);
    t = t ./ del;
    n = [-t(2) t(1)];
    correct_shape = size(xm);
    r = [(reshape(xm, [], 1) - xa) (reshape(ym, [], 1) - ya)];
    xm_rot = r * t';
    ym_rot = r * n';
    [infa, infb] = refpaninf(del, xm_rot, ym_rot);
    infa = reshape(infa, correct_shape);
    infb = reshape(infb, correct_shape);
end