clear;
clc;
close all;

%%          UNCERTAINTY MODEL TEST
 
% Description: test the uncetainty model derived using only Urr tensor 
% component. 


% IMPORTS
addpath("../partials_eq/functions/");
addpath("//Users/sergiocollibars/Desktop/CSML/codes/QGG/data_files/")
addpath("../../QGG_gravEstim/data_files/")
addpath("../matlab_codes/partials_eq/functions/")

% INPUT
N = 77760;                            % number of points

% Gravity field
n_max = 6;
Nc = 0;
for j = 1:n_max+1
    Nc  = Nc + j;
end
Nc = Nc - 2;
Ns = Nc - n_max;
Nx = Nc + Ns;
path = "HARMCOEFS_BENNU_CD_1.txt";
[C, S, R, normalized] = readCoeff(path);

% Planet parameters
deltaGM_GM = 0;
GM = 5.1355;                                                % [m^3/s^2]
GM_hat = GM + deltaGM_GM * GM;
Trot = 4.297461*3600;                                       % rotational period [s]
Wt = 2*pi/Trot;                                             % ang valocity planet [rad/s]

% orbit radius
Nrad = 3;
radius = linspace(500,1100,Nrad);
tMin = 0;
tMax = 777600;
TIME = linspace(tMin, tMax, N);
pos = zeros(3, N);
lon = zeros(1, N);
lat = zeros(1, N);

e = 0;                       % eccentrycity
a = radius(1);               % semi major axis [m]
i = deg2rad(90);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(0);              % true anomaly [rad]


% Noise value
W = (6.32E-13)^2;                                            % noise variance
Wmat_inv = eye(8) * (1/W);

% uncertainty
sigma_numerical = zeros(Nc+Ns, length(radius));
sigma_analytical = zeros(Nc+Ns, length(radius));
sigma_RSS = zeros(n_max, length(radius));

sigma_analytical_RSS_0 = zeros(n_max, length(radius));
sigma_analytical_RSS_1 = zeros(n_max, length(radius));
sigma_analytical_RSS_2 = zeros(n_max, length(radius));
sigma_analytical_RSS_3 = zeros(n_max, length(radius));
sigma_analytical_RSS_4 = zeros(n_max, length(radius));

for r=1:length(radius)
    clc;
    disp('radius iteration: ' + string(r));
    % Information matrix values
    Ax_numerical = 0;
    Ax_analytical = 0;

    % initial conditions. ACI (polar orbit)
    a = radius(r);
    rho = a * (1 - e^2);         % orbital param [m]
    [r0, v0] = compute_initialCond(rho, e, i, f, omega, Omega, GM);
    X0 = [r0; v0];
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    [~, state] = ode113(@(t, x) EoM(t, x, GM), TIME, X0, options);
    state = state';

    for k =1:N % longitude 
        % ACI coordinates 
        r_ACI = state(1:3, k);
        omega = Wt * TIME(k);
        ACAF_ACI = rotationMatrix(0, 0, omega, [3, 1, 3]);
        r_ACAF = ACAF_ACI * r_ACI;
        pos(:, k) = r_ACAF;

        % altitude
        h = radius(r);  

        % map coordinates
        r_xy = sqrt(r_ACAF(1)^2 + r_ACAF(2)^2);
        phi = atan2(r_ACAF(3), r_xy);
        lambda = atan2(r_ACAF(2), r_ACAF(1));

        lon(k) = lambda;
        lat(k) = phi;

        % transform to ENU coords
        ACAF_ENU = [-sin(lambda), -sin(phi)*cos(lambda), cos(phi)*cos(lambda);...
                    cos(lambda), -sin(phi)*sin(lambda), cos(phi)*sin(lambda);...
                    0, cos(phi), sin(phi)];

        ENU_ACAF = ACAF_ENU';
        
        % numerical H value. Visibility matrix
        [H_numerical] = potentialGradient_Cnm(n_max, r_ACAF, R, GM_hat, ...
            ENU_ACAF, normalized);
% %         H_numerical = [H_numerical(9, :); H_numerical(8, :);...
% %             H_numerical(7,:); H_numerical(5,:); H_numerical(4,:)];              %components: rr, rphi, rlambda, phiphi, philambda (all independent)
% %         H_numerical = [H_numerical(2:9, 1:2),H_numerical(2:9, 5),...
% %             H_numerical(2:9, 9), H_numerical(2:9, 14), H_numerical(2:9, 20)] ;      %components: only zonal
            H_numerical = H_numerical(2:9, :);

        Ax_numerical = Ax_numerical + (H_numerical' * Wmat_inv * H_numerical);
    end
    val = diag(Ax_numerical);
    %Ax_numerical = diag(val);

    % compute covariance. Numerical
    P_numerical = inv(Ax_numerical);
    sigma_numerical(:, r) = sqrt(diag(P_numerical))';
    sigma_RSS(:, r) = computeRSS_coefErr(n_max, Nc, Ns, sigma_numerical(:, r), ...
        zeros(n_max+1), zeros(n_max+1));

    % compute covariance. Analytical
    n = 0;
    pol= (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
    sigma_analytical(1, r) = (h^3/GM) * 1/sqrt(N) * sqrt(W)/sqrt(pol);
    n = 2;
    m = 0;
    for j =2:Nc
        pol= (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
        sigma_analytical(j, r) = ((h^3)/GM) * ((h/R)^n) * 1/sqrt(N) *...
            sqrt(W)/sqrt(pol);
        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
    end

    % compute RSS covariance models
    n = 0;
    pol_0 = (n+1)^2 *(n+2)^2 ;
    pol_1 = n*(n+1)*(n+2)^2;
    pol_2 = ((n-1)*(n+1)*n*(n+2));
    pol_3= (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
    pol_4 = pol_1 + pol_2 + pol_0;
    dev = (pol_4 - pol_3)/2;
    %pol_3 = pol_3 + dev;

    sigma_analytical_RSS_0(1, r) = (h^3/GM) * 1/sqrt(N) * sqrt(W)/sqrt(pol_0) * sqrt(2*n + 1);
    sigma_analytical_RSS_1(1, r) = (h^3/GM) * 1/sqrt(N) * sqrt(W)/sqrt(pol_1) * sqrt(2*n + 1);
    sigma_analytical_RSS_2(1, r) = (h^3/GM) * 1/sqrt(N) * sqrt(W)/sqrt(pol_2) * sqrt(2*n + 1);
    sigma_analytical_RSS_3(1, r) = (h^3/GM) * 1/sqrt(N) * sqrt(W)/sqrt(pol_3) * sqrt(2*n + 1);
    sigma_analytical_RSS_4(1, r) = (h^3/GM) * 1/sqrt(N) * sqrt(W)/sqrt(pol_4) * sqrt(2*n + 1);
    for j =2:n_max
        n = j;
        pol_0 = (n+1)^2 *(n+2)^2 ;
        pol_1 = n*(n+1)*(n+2)^2;
        pol_2 = ((n-1)*(n+1)*n*(n+2));
        pol_3= (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
        pol_4 = pol_1 + pol_2 + pol_0;
        dev = (pol_4 - pol_3)/2;
        %pol_3 = pol_3 + dev;
    
        sigma_analytical_RSS_0(j, r) = ((h^3)/GM) * ((h/R)^n) * 1/sqrt(N) *...
            sqrt(W)/sqrt(pol_0) * sqrt(2*n + 1);
        sigma_analytical_RSS_1(j, r) = ((h^3)/GM) * ((h/R)^n) * 1/sqrt(N) *...
            sqrt(W)/sqrt(pol_1) * sqrt(2*n + 1);
        sigma_analytical_RSS_2(j, r) = ((h^3)/GM) * ((h/R)^n) * 1/sqrt(N) *...
            sqrt(W)/sqrt(pol_2) * sqrt(2*n + 1);
        sigma_analytical_RSS_3(j, r) = ((h^3)/GM) * ((h/R)^n) * 1/sqrt(N) *...
            sqrt(W)/sqrt(pol_3) * sqrt(2*n + 1);
        sigma_analytical_RSS_4(j, r) = ((h^3)/GM) * ((h/R)^n) * 1/sqrt(N) *...
            sqrt(W)/sqrt(pol_4) * sqrt(2*n + 1);
        
    end
end


% PLOTS

% plot 3d ACI position
figure()
plot3(state(1, :), state(2, :), state(3, :), '*');
axis equal;
grid on;
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('ACI position')

%plot 3D ACAF position
figure()
plot3(pos(1, :), pos(2, :), pos(3, :), '*');
axis equal;
grid on;
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('ACAF position')

% plot ground track
figure()
plot(rad2deg(lon), rad2deg(lat), 'Marker','diamond','MarkerFaceColor', ...
    'magenta', 'LineStyle','none');
ylabel('Latitude, \phi [deg]')
xlabel('Longitude, \lambda [deg]')
title('Ground track sampling')

% plot RSS sigma vs value anlytical
[p, ~] = numSubplots(n_max);
rows = p(1);
cols = p(2);
figure();
subplot(rows, cols, 1)
plot(radius, sigma_RSS(1, :), '*', 'Color','magenta')
hold all;
plot(radius, sigma_analytical_RSS_0(1, :), ...
        'LineWidth', 1.5);
plot(radius, sigma_analytical_RSS_1(1, :), ...
        'LineWidth', 1.5);
plot(radius, sigma_analytical_RSS_2(1, :), ...
        'LineWidth', 1.5);
plot(radius, sigma_analytical_RSS_3(1, :), ...
    'LineWidth', 1.5);
plot(radius, sigma_analytical_RSS_4(1, :), ...
    'LineWidth', 1.5);
title('C_{00}');
% % legend('Numerical', '\alpha = 0', '\alpha = 1', '\alpha = 2', ...
% %     'UPM', '\alpha = \alpha_0 + \alpha_1 + \alpha_2')
legend('Numerical', '\alpha = 0', '\alpha = 1', '\alpha = 2', ...
    'UPM', '\alpha 4')
for j = 2:n_max
    subplot(rows, cols, j)
    plot(radius, sigma_RSS(j, :), '*', 'Color','magenta')
    hold all;
    plot(radius, sigma_analytical_RSS_0(j, :), ...
         'LineWidth', 1.5);
    plot(radius, sigma_analytical_RSS_1(j, :), ...
        'LineWidth', 1.5);
    plot(radius, sigma_analytical_RSS_2(j, :), ...
         'LineWidth', 1.5);
    plot(radius, sigma_analytical_RSS_3(j, :), ...
        'LineWidth', 1.5);
    plot(radius, sigma_analytical_RSS_4(j, :), ...
        'LineWidth', 1.5);
end 
sgtitle('RSS \sigma vs analytical model')

% plot RSS sigma vs degree
figure();
names = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", ...
    "14", "15", "16"];
lw = 1.5;
semilogy(linspace(2, n_max, n_max-1), sigma_RSS(2:end, :), 'sq--', ...
    'LineWidth', lw, 'Color', 'k')
hold all;
% % semilogy(linspace(2, n_max, n_max-1), sigma_analytical_RSS_0(2:end, :), '-', ...
% %     'LineWidth', lw)
% % semilogy(linspace(2, n_max, n_max-1), sigma_analytical_RSS_1(2:end, :), '-', ...
% %     'LineWidth', lw)
% % semilogy(linspace(2, n_max, n_max-1), sigma_analytical_RSS_2(2:end, :), '-', ...
% %     'LineWidth', lw)
loglog(linspace(2, n_max, n_max-1), sigma_analytical_RSS_3(2:end, :), '-', ...
    'LineWidth', lw)
loglog(linspace(2, n_max, n_max-1), sigma_analytical_RSS_4(2:end, :), '-', ...
    'LineWidth', lw)
%legend('Numerical', '\alpha = 0', '\alpha = 1', '\alpha = 2', '\alpha = 3')
legend('Numerical', '\alpha = 3', '\alpha = 4')
ylabel('\sigma RSS');
xlabel('Degree, n');
set(gca, 'XTick',2:n_max, 'XTickLabel',names)
title('Numerical vs Analytical degree uncertainty');

% plot individual sigmas vs value analytical
[p, ~] = numSubplots(Nc);
rows = p(1);
cols = p(2);
figure()
subplot(rows, cols, 1)
plot(radius, sigma_numerical(1, :), '*', 'Color','magenta')
hold on;
plot(radius, sigma_analytical(1, :), 'k', 'LineWidth', 1.5);
xlabel('radius')
ylabel('\sigma')
title("C_{00}");

degree = 2;
order = 0;
for k =2:Nc
    subplot(rows, cols, k)
    plot(radius, sigma_numerical(k, :), '*', 'Color','magenta')
    hold on;
    plot(radius, sigma_analytical(k, :), 'k', 'LineWidth', 1.5);
    xlabel('radius')
    ylabel('\sigma')
    name = "C_{" + string(degree) + string(order) + '}';
    title(name);
    if(order < degree)
        order = order + 1;
    else
        degree = degree + 1;
        order = 0;
    end
    title(string(name))
end


%% FUNCTION
function [CoefErr] = computeRSS_coefErr(n_max, Nc, Ns, X, Cnm, Snm)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = (X(1))^2;
    m  = 0;
    n = 2;
    for j =2:Nc
        CoefErr(n) = CoefErr(n)  + (X(j) - Cnm(n + 1, m+1))^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        CoefErr(n) = CoefErr(n)  + (X(j) - Snm(n + 1, m+1))^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
    for n = 1:n_max
        CoefErr(n) = sqrt(CoefErr(n));
    end
end

function [corr] = compute_correlation(Nc, Ns, P_numerical)
    corr = zeros(Nc+Ns, Nc + Ns);
    for i = 1:Ns+Nc
        for j = 1:Ns+Nc
            corr(i, j) = abs(P_numerical(i, j) / (sqrt(P_numerical(i, i))...
                *sqrt(P_numerical(j, j))));
        end
    end
end
function [dx] = EoM(t, x, GM)
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    dU = GM * [-x(1)/r^3, -x(2)/r^3, -x(3)/r^3];
    dx = [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3)];
    
end

function [p,n]=numSubplots(n)
    % function [p,n]=numSubplots(n)
    %
    % Purpose
    % Calculate how many rows and columns of sub-plots are needed to
    % neatly display n subplots. 
    %
    % Inputs
    % n - the desired number of subplots.     
    %  
    % Outputs
    % p - a vector length 2 defining the number of rows and number of
    %     columns required to show n plots.     
    % [ n - the current number of subplots. This output is used only by
    %       this function for a recursive call.]
    %
    %
    %
    % Example: neatly lay out 13 sub-plots
    % >> p=numSubplots(13)
    % p = 
    %     3   5
    % for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
    %
    %
    % Rob Campbell - January 2010
       
        
    while isprime(n) & n>4, 
        n=n+1;
    end
    p=factor(n);
    if length(p)==1
        p=[1,p];
        return
    end
    while length(p)>2
        if length(p)>=4
            p(1)=p(1)*p(end-1);
            p(2)=p(2)*p(end);
            p(end-1:end)=[];
        else
            p(1)=p(1)*p(2);
            p(2)=[];
        end    
        p=sort(p);
    end
    %Reformat if the column/row ratio is too large: we want a roughly
    %square design 
    while p(2)/p(1)>2.5
        N=n+1;
        [p,n]=numSubplots(N); %Recursive!
    end
end

function [r0, v0] = compute_initialCond(rho, e, i, f, omega, Omega, GM)
    r0 = rho / (1 + e * cos(f)) * [cos(f);...
                                      sin(f);...
                                      0];
    v0 = sqrt(GM / rho) * [-sin(f);...
                          e + cos(f);...
                          0];
    [BN] = rotationMatrix(Omega, i, omega, [3,1,3]);
    r0 = BN' * r0;
    v0 = BN' * v0;
end
  