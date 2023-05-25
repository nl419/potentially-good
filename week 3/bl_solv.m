function [int, ils, itr, its, delstar, theta] = bl_solv(x,cp)
    global Re ue0
    n_trapeziums = 100;
    ue = sqrt(1 - cp);

    theta = zeros([n_trapeziums + 1,1]);
    He = zeros([n_trapeziums + 1,1]);
    He(1) = 1.57258;

    laminar = true;
    i = 1;

    int = 0;
    ils = 0;
    itr = 0;
    its = 0;

    result = ueintbit(0, 0, x(1), ue(1));

    while laminar && i < (n_trapeziums+1)
        i = i + 1;
        result = result + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        theta(i) = sqrt((0.45 / Re) * power(ue(i), -6) * result);
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

    while its == 0 && i < (n_trapeziums+1)
        i = i + 1;
        thick0 = [theta(i-1), de(i-1)];
        ue0 = ue(i);
        [delx, thickhist] = ode45(@thickdash, [0, x(i) - x(i-1)], thick0);
        theta(i) = thickhist(end,1);
        de(i) = thickhist(end,2);
        He(i) = thickhist(end,2) / thickhist(end,1);

        if He(i) >= 1.46
            H = (11 * He + 15) / (48 * He - 59);
        else
            H = 2.803;
        end
        delstar(i) = H * theta(i);

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
        He(i+1) = He(i);
        theta(i) = theta(i-1) * power(ue(i-1) / ue(i), H + 2);
        delstar(i) = H * theta(i);
    end
end