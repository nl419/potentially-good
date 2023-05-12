function lhsmat = build_lhs(xs,ys)
    np = length(xs) - 1;
    psip = zeros(np,np+1); 
    for i = 1:np
        for  j = 1:(np+1)
           [fa1, fb1] = panelinf(xs(j), ys(j), ...
        xs(j+1), ys(j+1), xs(i), ys(i));
            if j == 1
                psi(i,j) = fa1;
            elseif j == (np+1)
                psi(i, j) = fb1;
            else
                psi(i, j) = fa1 + fb1;
            end
        end
    end
    
    lhsmat = zeros(np+1,np+1);
    
end

