function lhsmat = build_lhs(xs,ys)
    np = length(xs) - 1;
    psip = zeros(np+1,np+1);
    % This loop is the bottleneck of our entire program
    for i = 1:(np+1)
        [fa1, fb1] = panelinf(xs(np), ys(np), ...
            xs(1), ys(1), xs(i), ys(i));
        psip(i, np+1) = fb1;
        [fa1, fb1] = panelinf(xs(1), ys(1), ...
            xs(2), ys(2), xs(i), ys(i));
        psip(i,1) = fa1;
        fb2 = fb1;
        for j = 2:np
            [fa1, fb2_next] = panelinf(xs(j), ys(j), ...
                xs(j+1), ys(j+1), xs(i), ys(i));
            psip(i, j) = fa1 + fb2;
            fb2 = fb2_next;
        end
    end
    lhsmat = zeros(np+1,np+1);
    psi_n2 = circshift(psip,-1,1);
    psi_t = psi_n2 - psip;
    lhsmat = psi_t;

    lhsmat(1,:) = zeros(1,np+1); % FIXME
    lhsmat(end,:) = zeros(1,np+1); % FIXME
    lhsmat(1,1) = 1;
    lhsmat(1,2) = -1;
    lhsmat(1,3) = 0.5;
    lhsmat(1,np) = 1;
    lhsmat(1,np-1) = -0.5;

    lhsmat(end, 1) = 1;
    lhsmat(end, end) = 1;
end

