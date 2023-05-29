function lhsmat = build_lhs(xs,ys)
    np = length(xs) - 1;
    psip = zeros(np+1,np+1);
    % This loop is the bottleneck of our entire program
    for i = 1:np
        [fa, fb] = panelinf(xs(i), ys(i), ...
            xs(i+1), ys(i+1), xs, ys);
        psip(:,i) = psip(:,i) + fa';
        psip(:,i+1) = psip(:,i+1) + fb';
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

