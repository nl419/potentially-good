function [xk yk] = nik10_nowrite(params)
    % geometry parameters from section id designator
    m = params(0);
    p = params(1);
    tmax = params(2);
    tc0
    

    % camber line and thickness
    xc = [0 .0015 .003 0.006 .012 .025:.025:.1 .15:.05:.3 .4:.1:.9 .95 1];
    nc = length(xc);
    yc(1) = 0;
    t(1) = 0;
    for ic = 2:nc-1
      if xc(ic)<=p 
        yc(ic) = m * (2*p - xc(ic)) * xc(ic)/p^2;
      else
        yc(ic) = m * (1 - 2*p + 2*p*xc(ic) - xc(ic)^2)/(1 - p)^2;
      end 
      t(ic) = tmax * (2.969*sqrt(xc(ic)) - 1.26*xc(ic) - 3.516*xc(ic)^2 ...
                                            + 2.843*xc(ic)^3 - 1.036*xc(ic)^4);
    end
    yc(nc) = 0;
    t(nc) = 0;
    
    % surface coordinates
    
    xk(1) = 1;
    yk(1) = 0;
    
    theta = atan ( (yc(3:nc) - yc(1:nc-2))./(xc(3:nc) - xc(1:nc-2)) );
    
    xk(nc-1:-1:2) = xc(2:nc-1) - 0.5*t(2:nc-1).*sin(theta);
    xk(nc+1:2*nc-2) = xc(2:nc-1) + 0.5*t(2:nc-1).*sin(theta);
    
    yk(nc-1:-1:2) = yc(2:nc-1) + 0.5*t(2:nc-1).*cos(theta);
    yk(nc+1:2*nc-2) = yc(2:nc-1) - 0.5*t(2:nc-1).*cos(theta);
    
    xk(2*nc-1) = xk(1);
    yk(2*nc-1) = yk(1);
    % Make the vector orientation consistent with the file we just wrote (Nik).
    xk = xk';
    yk = yk';
end
    