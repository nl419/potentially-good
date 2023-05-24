clc
close all
clear all

global Re ue0 duedx

%% Initial
figure;
Re_vals = [1e4 1e5 1e6];
duedx = -0.25;

%% Loops
for i = 1:length(Re_vals)
    Re = Re_vals(i);
    n_trapeziums = 100;
    x = linspace(0,1,n_trapeziums + 1);
    ue = ones([n_trapeziums + 1,1]);
    for i=1:n_trapeziums + 1
        ue(i) = 1 + duedx * (i-1) / n_trapeziums;
    end

    result = 0;

    theta = zeros([n_trapeziums + 1,1]);
    He = zeros([n_trapeziums + 1,1]);

    laminar = true;
    i = 1;

    int = 0;
    ils = 0;
    itr = 0;
    its = 0;
    Rethet = 0;

    while laminar && i < (n_trapeziums+1)
        i = i + 1;
        result = result + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        theta(i) = sqrt((0.45 / Re) * power(ue(i), -6) * result);
        m = -Re * theta(i) * theta(i) * duedx;
        H = thwaites_lookup(m);
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
    end

    if int ~= 0
        i = int;
        disp(['Natural transition at ' num2str(x(i)) ' with Rethet ' num2str(Re * ue(i) * theta(i))]);
    end
    if ils ~= 0
        i = ils;
        disp(['Laminar separation at ' num2str(x(i)) ' with Rethet ' num2str(Re * ue(i) * theta(i))]);
    end
    if itr ~= 0
        i = itr;
        disp(['Turbulent reattachment at ' num2str(x(i)) ' with Rethet ' num2str(Re * ue(i) * theta(i))]);
    end
    if its ~= 0
        i = its;
        disp(['Turbulent separation at ' num2str(x(i)) ' with Rethet ' num2str(Re * ue(i) * theta(i))]);
    end
    scatter(x, He);
    hold on
end

hold off
legend({'$Re_L = 1 \times 10^6$', '$Re_L = 1 \times 10^7$'}, 'interpreter','latex','FontSize',16)
xlabel("x/L", 'FontSize',18,'FontWeight','bold');
ylabel("$H_E$", 'interpreter','latex', 'FontSize',18,'FontWeight','bold');
% xline(0.37,'--r',{'Natural Transition'}, 'HandleVisibility','off')
xline(0.5,'--b',{'Laminar Separation'}, 'HandleVisibility','off')
xline(0.5,'--r',{'Laminar Separation'}, 'HandleVisibility','off')
xline(0.59,'--r',{'Turbulent Reattachment'}, 'HandleVisibility','off')
xline(0.49,'--g',{'Natural Transition'}, 'HandleVisibility','off')
ylim([1.45 1.9]);