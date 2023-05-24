clc
close all
clear all

global Re ue0 duedx

%% Initial
figure;
Re_vals = [1e4 1e5 1e6];
gradU_norm = -0.25;

%% Loops
for i = 1:length(Re_vals)
    Re = Re_vals(i);
    n_trapeziums = 100;
    xs_norm = linspace(0,1,n_trapeziums + 1);
    us_norm = ones([n_trapeziums + 1,1]);
    for i=1:n_trapeziums + 1
        us_norm(i) = 1 + gradU_norm * (i-1) / n_trapeziums;
    end

    result = 0;

    theta_over_L = zeros([n_trapeziums + 1,1]);
    He = zeros([n_trapeziums + 1,1]);

    laminar = true;
    i = 0;

    int = 0;
    ils = 0;
    itr = 0;
    its = 0;
    Rethet = 0;

    while laminar && i < n_trapeziums
        i = i + 1;
        result = result + ueintbit(xs_norm(i), us_norm(i), xs_norm(i+1), us_norm(i+1));
        theta_over_L(i+1) = sqrt((0.45 / Re) * power(us_norm(i+1), -6) * result);
        m = -Re * theta_over_L(i+1) * theta_over_L(i+1) * gradU_norm;
        H = thwaites_lookup(m);
        He(i+1) = laminar_He(H);
        Rethet = Re * us_norm(i+1) * theta_over_L(i+1);

        if log(Rethet) >= 18.4*He(i+1) - 21.74
            laminar = false;
            int = i+1;
        elseif m >= 0.09
            laminar = false;
            ils = i+1;
            He(i+1) = 1.51509;
        end
    end

    de_over_L = He .* theta_over_L;

    while its == 0 && i < n_trapeziums
        i = i + 1;
        thick0 = [theta_over_L(i), de_over_L(i)];
        ue0 = us_norm(i);
        duedx = gradU_norm;
        [delx, thickhist] = ode45(@thickdash, [0, xs_norm(i+1) - xs_norm(i)], thick0);
        theta_over_L(i+1) = thickhist(end,1);
        de_over_L(i+1) = thickhist(end,2);
        He(i+1) = thickhist(end,2) / thickhist(end,1);

        % Test for turb sep always
        if its == 0 && itr ~= 0 && He(i+1) < 1.46
            its = i+1;
        end
        % Test for reattachment if ils != 0
        if ils ~= 0 && itr == 0 && He(i+1) > 1.58
            itr = i+1;
        end
    end

    while its ~= 0 && i < n_trapeziums
        i = i + 1;
        H = 2.803;
        He(i+1) = He(i);
        theta_over_L(i+1) = theta_over_L(i) * power(us_norm(i) / us_norm(i + 1), H + 2);
    end

    if int ~= 0
        i = int;
        disp(['Natural transition at ' num2str(xs_norm(i)) ' with Rethet ' num2str(Re * us_norm(i) * theta_over_L(i))]);
    end
    if ils ~= 0
        i = ils;
        disp(['Laminar separation at ' num2str(xs_norm(i)) ' with Rethet ' num2str(Re * us_norm(i) * theta_over_L(i))]);
    end
    if itr ~= 0
        i = itr;
        disp(['Turbulent reattachment at ' num2str(xs_norm(i)) ' with Rethet ' num2str(Re * us_norm(i) * theta_over_L(i))]);
    end
    if its ~= 0
        i = its;
        disp(['Turbulent separation at ' num2str(xs_norm(i)) ' with Rethet ' num2str(Re * us_norm(i) * theta_over_L(i))]);
    end
    scatter(xs_norm, He);
    hold on
end

hold off
legend({'$Re_L = 1 \times 10^4$', '$Re_L = 1 \times 10^5$', '$Re_L = 1 \times 10^6$'}, 'interpreter','latex','FontSize',16)
xlabel("x/L", 'FontSize',18,'FontWeight','bold');
ylabel("$\theta/L$", 'interpreter','latex', 'FontSize',18,'FontWeight','bold');
xline(0.5,'--b',{'Laminar Separation'}, 'HandleVisibility','off')
xline(0.5,'--r',{'Laminar Separation'}, 'HandleVisibility','off')
xline(0.59,'--r',{'Turbulent Reattachment'}, 'HandleVisibility','off')
xline(0.49,'--g',{'Natural Transition'}, 'HandleVisibility','off')
ylim([1 2]);