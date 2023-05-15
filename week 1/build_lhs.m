function lhsmat = build_lhs(xs,ys)
    np = length(xs) - 1;
    psip = zeros(np+1,np+1); 
    for i = 1:(np+1)
        for  j = 1:(np+1)
            if j == 1
                [fa1, fb1] = panelinf(xs(j), ys(j), ...
        xs(j+1), ys(j+1), xs(i), ys(i));
                psip(i,j) = fa1;
            elseif j == np+1
                [fa1, fb1] = panelinf(xs(j-1), ys(j-1), ...
        xs(1), ys(1), xs(i), ys(i));
                psip(i, j) = fb1;
            else
               [fa1, fb1] = panelinf(xs(j), ys(j), ...
        xs(j+1), ys(j+1), xs(i), ys(i));
               [fa2, fb2] = panelinf(xs(j-1), ys(j-1), ...
        xs(j), ys(j), xs(i), ys(i));
                psip(i, j) = fa1 + fb2;
            end
        end
    end
    
    lhsmat = zeros(np+1,np+1);
    psi_n2 = circshift(psip,-1,1);
    psi_t = psi_n2 - psip;
    lhsmat = psi_t; 
    lhsmat(end,:) = [];
end

