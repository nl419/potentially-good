function [int, ils, itr, its, delstar, theta] = bl_solv(x,cp,Re)
    n_trapeziums = length(x);
    ue = sqrt(1 - cp);

    theta = zeros([n_trapeziums,1]);
    He = zeros([n_trapeziums,1]);
    He(1) = 1.57258;

    laminar = true;
    i = 1;

    int = 0;
    ils = 0;
    itr = 0;
    its = 0;

    result = ueintbit(0, 0, x(1), ue(1));
    theta(1) = sqrt((0.45 / Re) * power(ue(1), -6) * result);

    while laminar && i < (n_trapeziums)
        i = i + 1;
        result = result + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        theta(i) = sqrt((0.45 / Re) * power(ue(i), -6) * result);
        duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1));
        m = -Re * theta(i) * theta(i) * duedx;
        H = thwaites_lookup(m);
        delstar(i) = H * theta(i);
        He(i) = laminar_He(H);
        Rethet = Re * ue(i) * theta(i);

        if log(Rethet) >= 18.4*He(i) - 21.74
            laminar = false;
            int = i;
        elseif m >= 0.09
            laminar = false;
            ils = i;
            He(i) = 1.51509;
        end
    end

    de = He .* theta;

    while its == 0 && i < (n_trapeziums)
        i = i + 1;
        thick0 = [theta(i-1), de(i-1)];
        duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1));
        ue0 = ue(i);
        [delx, thickhist] = ode45(@(x,y) thickdash(x,y,Re,ue0,duedx), [0, x(i) - x(i-1)], thick0);
        theta(i) = thickhist(end,1);
        de(i) = thickhist(end,2);
        He(i) = thickhist(end,2) / thickhist(end,1);

        if He(i) >= 1.46
            H = (11 * He(i) + 15) / (48 * He(i) - 59);
            delstar(i) = H * theta(i);
        else
            H = 2.803;
            delstar(i) = H * theta(i);
        end


        % Test for turb sep always
        if its == 0 && itr ~= 0 && He(i) < 1.46
            its = i;
        end
        % Test for reattachment if ils != 0
        if ils ~= 0 && itr == 0 && He(i) > 1.58
            itr = i;
        end
    end

    while its ~= 0 && i < n_trapeziums
        i = i + 1;
        H = 2.803;
        He(i) = He(i-1);
        theta(i) = theta(i-1) * power(ue(i-1) / ue(i), H + 2);
        delstar(i) = H * theta(i);
    end
end