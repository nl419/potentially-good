function [rhsvec] = build_rhs(xs,ys,alpha)
    np = length(xs) - 1;
    b = zeros(np+1,1);
    psi_n = ys.*cos(alpha) - xs.*sin(alpha);
    psi_n2 = circshift(psi_n,-1);
    psi_t = psi_n - psi_n2;
    rhsvec = psi_t';
    rhsvec(1) = 0;
    rhsvec(end) = 0;
end