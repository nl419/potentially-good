close all
clear all
clc

global Re

np = 400;
Re = 20e6;
alpha = 0:0.2:10;

folderPath = 'Barty/'; % Include trailing `/`
fileExtension = '.surf';

% Get a list of all files in the folder
fileList = dir(folderPath);

% Remove directories from the file list
fileList = fileList(~[fileList.isdir]);

% Initialize variables
mostRecentFile = '';
mostRecentDate = 0;

% Iterate through the files and find the most recently modified file with the specified extension
for i = 1:numel(fileList)
    [~, ~, extension] = fileparts(fileList(i).name);
    
    % Check if the file has the desired extension
    if strcmpi(extension, fileExtension)
        % Check if the file's modification date is more recent
        if fileList(i).datenum > mostRecentDate
            mostRecentDate = fileList(i).datenum;
            mostRecentFile = fileList(i).name;
        end
    end
end

if isempty(mostRecentFile)
    disp(['No files with extension ' fileExtension ' found in the folder.']);
    return;
else
    % Open the most recently modified file
    secfile = fullfile(folderPath, mostRecentFile);

    % Display the path of the most recently modified file
    disp(['Most recently modified file with extension ' fileExtension ': ' mostRecentFile]);
end

[xk yk] = textread ( secfile, '%f%f' );
[~, caseref, ~] = fileparts(secfile);

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshr yshr] = splinefit ( xk, yk, nphr );

%  Resize section so that it lies between (0,0) and (1,0)
[xsin ysin] = resyze ( xshr, yshr );

%  Interpolate to required number of panels (uniform size)
[xs ys] = make_upanels ( xsin, ysin, np );

%  Assemble the lhs of the equations for the potential flow calculation
A = build_lhs ( xs, ys );
Am1 = inv(A);

cpu_max = [];
cpl_max = [];
L_over_D_max = 0;
su_max = [];
sl_max = [];
ilnt_max = 0;
ills_max = 0;
iltr_max = 0;
ilts_max = 0;
iunt_max = 0;
iuls_max = 0;
iutr_max = 0;
iuts_max = 0;
alpha_max = 0;

%  Loop over alpha values
for nalpha = 1:length(alpha)
    %    rhs of equations
    alfrad = pi * alpha(nalpha)/180;
    b = build_rhs ( xs, ys, alfrad );

    %    solve for surface vortex sheet strength
    gam = Am1 * b;

    %    calculate cp distribution and overall circulation
    [cp circ] = potential_op ( xs, ys, gam );

    %    locate stagnation point and calculate stagnation panel length
    [ipstag fracstag] = find_stag(gam);
    dsstag = sqrt((xs(ipstag+1)-xs(ipstag))^2 + (ys(ipstag+1)-ys(ipstag))^2);

    %    upper surface boundary layer calc

    %    first assemble pressure distribution along bl
    clear su cpu
    su(1) = fracstag*dsstag;
    cpu(1) = cp(ipstag);
    for is = ipstag-1:-1:1
        iu = ipstag - is + 1;
        su(iu) = su(iu-1) + sqrt((xs(is+1)-xs(is))^2 + (ys(is+1)-ys(is))^2);
        cpu(iu) = cp(is);
    end

    %    check for stagnation point at end of stagnation panel
    if fracstag < 1e-6
        su(1) = 0.01*su(2);    % go just downstream of stagnation
        uejds = 0.01 * sqrt(1-cpu(2));
        cpu(1) = 1 - uejds^2;
    end

    %    boundary layer solver
    [iunt iuls iutr iuts delstaru thetau] = bl_solv(su, cpu, Re);

    %    lower surface boundary layer calc

    %    first assemble pressure distribution along bl
    clear sl cpl
    sl(1) = (1-fracstag) * dsstag;
    cpl(1) = cp(ipstag+1);
    for is = ipstag+2:np+1
        il = is - ipstag;
        sl(il) = sl(il-1) + sqrt((xs(is-1)-xs(is))^2 + (ys(is-1)-ys(is))^2);
        cpl(il) = cp(is);
    end

    %    check for stagnation point at end of stagnation panel
    if fracstag > 0.999999
        sl(1) = 0.01*sl(2);    % go just downstream of stagnation
        uejds = 0.01 * sqrt(1-cpl(2));
        cpl(1) = 1 - uejds^2;
    end

    %    boundary layer solver
    [ilnt ills iltr ilts delstarl thetal] = bl_solv(sl, cpl, Re);

    %    lift and drag coefficients
    [Cl Cd] = forces ( circ, cp, delstarl, thetal, delstaru, thetau );

    %    copy Cl and Cd into arrays for alpha sweep plots

    clswp(nalpha) = Cl;
    cdswp(nalpha) = Cd;
    lovdswp(nalpha) = Cl/Cd;

    L_over_D = Cl/Cd;

    if L_over_D > L_over_D_max 
        if length(alpha) ~= 1 && iuts ~= 0 && su(iuts) < 0.9
        else
            L_over_D_max = L_over_D;
            cpu_max = cpu;
            cpl_max = cpl;
            su_max = su;
            sl_max = sl;
            iunt_max = iunt;
            iuls_max = iuls;
            iutr_max = iutr;
            iuts_max = iuts;
            ilnt_max = ilnt;
            ills_max = ills;
            iltr_max = iltr;
            ilts_max = ilts;
            alpha_max = alpha(nalpha);
        end
    end

    %    screen output

    disp ( sprintf ( '\n%s%5.3f%s', ...
                    'Results for alpha = ', alpha(nalpha), ' degrees' ) )

    disp ( sprintf ( '\n%s%5.3f', '  Lift coefficient: ', Cl ) )
    disp ( sprintf ( '%s%7.5f', '  Drag coefficient: ', Cd ) )
    disp ( sprintf ( '%s%5.3f\n', '  Lift-to-drag ratio: ', Cl/Cd ) )

    upperbl = sprintf ( '%s', '  Upper surface boundary layer:' );
    if iunt~=0
        is = ipstag + 1 - iunt;
        upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                            '    Natural transition at x = ', xs(is) );
    end
    if iuls~=0
        is = ipstag + 1 - iuls;
        upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                            '    Laminar separation at x = ', xs(is) );
        if iutr~=0
        is = ipstag + 1 - iutr;
        upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                            '    Turbulent reattachment at x = ', xs(is) );
        end
    end
    if iuts~=0
        is = ipstag + 1 - iuts;
        upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                            '    Turbulent separation at x = ', xs(is) );
    end
    upperbl = sprintf ( '%s\n', upperbl );
    disp(upperbl)

    lowerbl = sprintf ( '%s', '  Lower surface boundary layer:' );
    if ilnt~=0
        is = ipstag + ilnt;
        lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                            '    Natural transition at x = ', xs(is) );
    end
    if ills~=0
        is = ipstag + ills;
        lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                            '    Laminar separation at x = ', xs(is) );
        if iltr~=0
        is = ipstag + iltr;
        lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                            '    Turbulent reattachment at x = ', xs(is) );
        end
    end
    if ilts~=0
        is = ipstag + ilts;
        lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                            '    Turbulent separation at x = ', xs(is) );
    end
    lowerbl = sprintf ( '%s\n', lowerbl );
    disp(lowerbl);

    %    save data for this alpha
    % fname = [folderPath caseref '_' num2str(alpha(nalpha)) '.mat'];
    % save ( fname, 'Cl', 'Cd', 'xs', 'cp', ...
            % 'sl', 'delstarl', 'thetal', 'lowerbl', ...
            % 'su', 'delstaru', 'thetau', 'upperbl' )
end

% scatter(clswp, cdswp);
figure;
scatter(alpha, lovdswp);
title('L/D vs alpha');
xlabel('alpha');
ylabel('L/D');

if iunt_max ~= 0
    xline(su_max(iunt_max),'--b',{'Natural Transition'}, 'HandleVisibility','off');
end
if iuls_max ~= 0
    xline(su_max(iuls_max),'--b',{'Laminar Separation'}, 'HandleVisibility','off');
end
if iutr_max ~= 0
    xline(su_max(iutr_max),'--b',{'Turbulent Reattachment'}, 'HandleVisibility','off');
end
if iuts_max ~= 0
    xline(su_max(iuts_max),'--r',{'Turbulent Separation'}, 'HandleVisibility','off');
end

figure;
plot(sl_max, -cpl_max);
title(['Pressure coefficient on lower surface at alpha=' num2str(alpha_max)]);
xlabel('X');
ylabel('-Cp');
if ilnt_max ~= 0
    xline(sl_max(ilnt_max),'--b',{'Natural Transition'}, 'HandleVisibility','off');
end
if ills_max ~= 0
    xline(sl_max(ills_max),'--b',{'Laminar Separation'}, 'HandleVisibility','off');
end
if iltr_max ~= 0
    xline(sl_max(iltr_max),'--b',{'Turbulent Reattachment'}, 'HandleVisibility','off');
end
if ilts_max ~= 0
    xline(sl_max(ilts_max),'--r',{'Turbulent Separation'}, 'HandleVisibility','off');
end

store_mat = [L_over_D_max, alpha_max];
disp(store_mat);
