clc;
close all;
clear;
format long g;

%%                          SDIM FUNCTION TEST
% Description: test script for the State Deviation Inclusion Method (SDIM)
%   algorithm.

% INPUT
addpath('functions/');
addpath('dataFiles/');
perturbation = "bias"; % rnd = random / bias / sinusoidal
coefErr_truth  = load("dataFiles/CoeffErr_truth.mat").CoefErr;
coefErr_truth_noise_E12  = load("dataFiles/CoeffErr_truth_noise_E12.mat").CoefErr;

% Harmonics values. Truth
Cnm_t = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.0391557863539988 -2.96928723209235e-06 0.00375640654748657 0 0 0 0;...
        0.0148427177700986 0.00167095097673949 3.80845003468165e-05 0.000371755938456641 0 0 0;...
        0.0307494000000000 0.000413625917950024 -0.000490123739988179 -6.43092753252371e-05 4.51227856599555e-05 0 0;...
        -0.00456599734888228 0.000441961635589938 8.84264903209225e-05 1.40396087936725e-05 1.07458402471302e-05 1.19886311472876e-06 0;...
        -0.00896736657720649 0.000905916675449317 9.05780703119059e-05 2.77025499573633e-05 -1.92680576794137e-06 9.97533032070693e-08 1.67034838314692e-07];

Snm_t = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 0 -2.54325906400954e-05 0 0 0 0;...
        0 0.000992000134408593 6.53000000000000e-05 -0.00100797120329237 0 0 0;
        0 0.000634013000392474 0.000108054642426876 0.000102400000000000 0.00291093983173820 0 0;...
        0 -1.75754943031483e-05 -7.41878397813858e-05 4.79413750994751e-06 -0.000503800000000000 0.000448812426298560 0;...
        0 -1.35941163743731e-06 -9.69889209840526e-06 -7.55736000569855e-06 -1.59676058718751e-06 -0.000599300000000000 -3.93397896234722e-05];

% Harmonic values. A priori
Cnm_ap = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.04 -3e-06 0.004 0 0 0 0;...
        0.014 0.00179 3-05 0.0001 0 0 0;...
        0.02 0.0001 -0.0004 -7-05 4.51-05 0 0;...
        -0.005 0.00044 9-05 1e-05 1e-05 1e-06 0;...
        -0.009 0.001 9e-05 2e-05 -2e-06 9e-08 2e-07];

Snm_ap = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 1E-6 -3e-05 0 0 0 0;...
        0 0.0001 7e-05 -0.001 0 0 0;
        0 0.007 0.0001 0.0001 0.003 0 0;...
        0 -2e-05 -8e-05 5e-06 -0.0005 0.0005 0;...
        0 -1e-06 -1e-05 -8e-06 -2e-06 -0.0006 -4e-05];

% Bennu parameters
GM = 5.2;
Re = 246;
n = sqrt(GM / 1000^3);
T = 2 * pi / n;
n_max = 6;
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% number of points
t0 = 0;
t_max = round(777600);
Nt_max = round(t_max / 5);
t = linspace(t0, t_max, Nt_max);

% Nominal orbit
X0 = [3.69978448523722;
    -105.053702653589;
    994.459667937083;
    -0.0528661190622108;
    -0.0487913496705586;
    -0.00495758536239498];
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, GM, Re, Cnm_t, Snm_t, n_max, 0, W, W0, ...
    RA, DEC), t, X0, options);
state = state';

% Position error. RTN
delta_bias= [1; 1; 1];              % [m]

mean_rnd = [1; 0; 0];     % [m]
sigma_rnd = [1; 0; 0];    % [m]

amp_sinu = [1; 1; 1];     % [m]
omega_sinu = n;           % [rad/s]
bias_sinu  = 1;           % [m]

% nominal orbit. ACI
r_ACI = zeros(3, Nt_max);
r_ACAF = zeros(3, Nt_max);
r_RTN = zeros(3, Nt_max);
rr_ACI = zeros(3, Nt_max);
rr_ACAF = zeros(3, Nt_max);
rr_RTN = zeros(3, Nt_max);

% measurements vector
Y_o = zeros(9, Nt_max);
Y_oc = zeros(9, Nt_max);

% rotation matrix
ACAF_RTN_mat = zeros(3*Nt_max, 3);
RTN_ACI_mat = zeros(3*Nt_max, 3);
ACAF_ACI_mat = zeros(3*Nt_max, 3);

% generate measurements. Truth
sigma = 1E-12;
for k = 1:Nt_max
    rt_ACI = [state(1, k);state(2, k);state(3, k)];
    vt_ACI = [state(4, k);state(5, k);state(6, k)];
    r_ACI(:, k) = rt_ACI;
    
    % ACAF postion. Truth
    Wt = W0 + W * t(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    rt_ACAF = ACAF_ACI * rt_ACI;
    r_ACAF(:, k) = rt_ACAF;

     % RTN 2 ACI
    [ACI_RTN] = RTN2ECI(rt_ACI, vt_ACI);
    RTN_ACI = ACI_RTN';
    r_RTN(:,  k) = RTN_ACI * r_ACI(:, k);

    % RTN 2 ACAF
    ACAF_RTN = ACAF_ACI * ACI_RTN;
    up = 3*k;
    down = up - 2;
    ACAF_RTN_mat(down:up, :) = ACAF_RTN;
    RTN_ACI_mat(down:up, :) = RTN_ACI;
    ACAF_ACI_mat(down:up, :) = ACAF_ACI;

    % compute truth meas. RTN frame
    [~, ~, ddU] = potentialGradient_nm(Cnm_t, Snm_t, n_max, ...
                                           rt_ACAF, Re, GM, 0);
    ddU = ddU + sigma * randn(size(ddU));
    Y_o(:, k) = reshape(ACAF_RTN'* ddU * ACAF_RTN, [9, 1]);
end
Y_o = load('Ydata.mat').Y_o;

% inital values
X = [1; Cnm_ap(3,1); Cnm_ap(3,2); Cnm_ap(3,3); ...
    Cnm_ap(4,1); Cnm_ap(4,2); Cnm_ap(4,3); Cnm_ap(4,4);...
    Cnm_ap(5,1); Cnm_ap(5,2); Cnm_ap(5,3); Cnm_ap(5,4); Cnm_ap(5,5);...
    Cnm_ap(6,1); Cnm_ap(6,2); Cnm_ap(6,3); Cnm_ap(6,4); Cnm_ap(6,5); Cnm_ap(6,6);...
    Cnm_ap(7,1); Cnm_ap(7,2); Cnm_ap(7,3); Cnm_ap(7,4); Cnm_ap(7,5); Cnm_ap(7,6); Cnm_ap(7,7);...
    Snm_ap(3,2); Snm_ap(3,3); ...
    Snm_ap(4,2); Snm_ap(4,3); Snm_ap(4,4);...
    Snm_ap(5,2); Snm_ap(5,3); Snm_ap(5,4); Snm_ap(5,5);...
    Snm_ap(6,2); Snm_ap(6,3); Snm_ap(6,4); Snm_ap(6,5); Snm_ap(6,6);...
    Snm_ap(7,2); Snm_ap(7,3); Snm_ap(7,4); Snm_ap(7,5); Snm_ap(7,6); Snm_ap(7,7); ...
    0; 0; 0];

count = 1;
maxCount = 4;
errj = zeros(3*maxCount, Nt_max);
CoefErr = zeros(maxCount, n_max);
Xx = 0;
while(count <= maxCount)
    disp('Iteration: ' + string(count));
    % Batch values
    P0 = eye(48) * 1E1;
    if(count == 1)
        X0 = zeros(48, 1);
    else
        X0 = Xx;
    end
    %X0 = zeros(45, 1);
    Axp = inv(P0);
    Nxp = P0\-X0;
    
    % time loop
    for k = 1:Nt_max
        % ACAF RTN
        up = 3*k;
        down = up - 2;
        ACAF_RTN = ACAF_RTN_mat(down:up, :);
        
        % RTN ACI
        RTN_ACI = RTN_ACI_mat(down:up, :);

        % ACAF ACI
        ACAF_ACI = ACAF_ACI_mat(down:up, :);

        % ACAF position. Perturbed
        if(perturbation == "bias")
            rp_ACAF = r_ACAF(:, k) + (ACAF_RTN * delta_bias);
        elseif(perturbation == "rnd")
            delta = sigma_rnd.*rand(size(mean_rnd)) + mean_rnd;
            rp_ACAF = r_ACAF(:, k) + (ACAF_RTN * delta);
        elseif(perturbation == "sinusoidal")
           delta = amp_sinu * sin(omega_sinu * t(k)) + bias_sinu;
           rp_ACAF = r_ACAF(:, k) + (ACAF_RTN * delta);
        end

        % Correct position
        if(count == 1)
            rr_ACI(:, k) = ACAF_ACI' * rp_ACAF;
            rr_ACAF(:, k) = rp_ACAF;
        else % correct position
            rr_ACAF(:, k) = rr_ACAF(:, k) + ACAF_RTN * [XNOT(end-2); XNOT(end-1); XNOT(end)];
            rr_ACI(:, k) = ACAF_ACI' * rr_ACAF(:, k);
            rp_ACAF = rr_ACAF(:, k);
        end
        rr_RTN(:,  k) = RTN_ACI * rr_ACI(:, k);
    
        % compute perturbed meas. RTN frame
        [~, ~, ddU] = potentialGradient_nm(Cnm_ap, Snm_ap, n_max, ...
                                               rp_ACAF, Re, GM, 0);
        Yc = reshape(ACAF_RTN'* ddU * ACAF_RTN, [9, 1]);
    
        % H perturbed. ENU
        [Hp] = potentialGradient_Cnm(n_max, rp_ACAF, Re, GM, ACAF_RTN');
    
        % new meas
        Y = Y_o(:, k) - Yc;
        Y_oc(:, k) = Y;
        Y = [Y(1); Y(2); Y(3); Y(5); Y(6)];
        Hp = [Hp(1,2:end); Hp(2, 2:end); Hp(3, 2:end); Hp(5,2:end); Hp(6,2:end)];

        % include position error state
        rp = vecnorm(rp_ACAF);
        Hp_AR = 3*GM/rp^4 * [-2; 0; 0; 1; 0];
        Hp_AT = 3*GM/rp^4 * [0; 1; 0; 0; 0];
        Hp_AN = 3*GM/rp^4 * [0; 0; 1; 0; 0];
        Hp = [Hp, Hp_AR, Hp_AT, Hp_AN];
        R = eye(5, 5) * sigma^2;
        
        Axp = Axp + (Hp' * inv(R) * Hp);
        Nxp = Nxp + (Hp' * inv(R) * Y);
    end
    
    % Harmonic correction
    XNOT = Axp\Nxp;
    
    Cnm_updt = [1 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            XNOT(1)+Cnm_ap(3,1) XNOT(2)+Cnm_ap(3,2) XNOT(3)+Cnm_ap(3,3) 0 0 0 0;...
            XNOT(4)+Cnm_ap(4,1) XNOT(5)+Cnm_ap(4,2) XNOT(6)+Cnm_ap(4,3) XNOT(7)+Cnm_ap(4,4) 0 0 0;...
            XNOT(8)+Cnm_ap(5,1) XNOT(9)+Cnm_ap(5,2) XNOT(10)+Cnm_ap(5,3) XNOT(11)+Cnm_ap(5,4) XNOT(12)+Cnm_ap(5,5) 0 0;...
            XNOT(13)+Cnm_ap(6,1) XNOT(14)+Cnm_ap(6,2) XNOT(15)+Cnm_ap(6,3) XNOT(16)+Cnm_ap(6,4) XNOT(17)+Cnm_ap(6,5) XNOT(18)+Cnm_ap(6,6) 0;...
            XNOT(19)+Cnm_ap(7,1) XNOT(20)+Cnm_ap(7,2) XNOT(21)+Cnm_ap(7,3) XNOT(22)+Cnm_ap(7,4) XNOT(23)+Cnm_ap(7,5) XNOT(24)+Cnm_ap(7,6) XNOT(25)+Cnm_ap(7,7)];
    
    Snm_updt = [0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            0 XNOT(26)+Snm_ap(3,2) XNOT(27)+Snm_ap(3,3) 0 0 0 0;...
            0 XNOT(28)+Snm_ap(4,2) XNOT(29)+Snm_ap(4,3) XNOT(30)+Snm_ap(4,4) 0 0 0;
            0 XNOT(31)+Snm_ap(5,2) XNOT(32)+Snm_ap(5,3) XNOT(33)+Snm_ap(5,4) XNOT(34)+Snm_ap(5,5) 0 0;...
            0 XNOT(35)+Snm_ap(6,2) XNOT(36)+Snm_ap(6,3) XNOT(37)+Snm_ap(6,4) XNOT(38)+Snm_ap(6,5) XNOT(39)+Snm_ap(6,6) 0;...
            0 XNOT(40)+Snm_ap(7,2) XNOT(41)+Snm_ap(7,3) XNOT(42)+Snm_ap(7,4) XNOT(43)+Snm_ap(7,5) XNOT(44)+Snm_ap(7,6) XNOT(45)+Snm_ap(7,7)];
    
    % update state
    X(2:end) = X(2:end) + XNOT;

    Cnm_ap = Cnm_updt;
    Snm_ap = Snm_updt;
    
    Xx = Xx + XNOT;
    CoefErr(count, :) = computeRMS_coefErr(n_max, 26, 20, X(1:end-3), ...
        Cnm_t, Snm_t);

    % compute error
    up = count*3;
    down = up - 2;
    errj(down:up, :) = r_RTN - rr_RTN;

    % update counter
    count = count + 1;
end

%% PLOTS
% covariance analysis
P = inv(Axp);
sigma2 = diag(P);
sigma = sqrt(sigma2);

corr = zeros(48, 48);
for i = 1:48
    for j = 1:48
        corr(i, j) = P(i, j) / (sqrt(P(i, i))*sqrt(P(j, j)));
    end
end
Axcorr = zeros(1, 45);
Aycorr = zeros(1, 45);
Azcorr = zeros(1, 45);
[rx,ry,rz] = deal(1,1,1);
for j = 1:48
    if(j ~= 46)
        Axcorr(rx) = corr(46, j);
        rx = rx + 1;
    end
    if(j ~= 47)
        Aycorr(ry) = corr(47, j);
        ry = ry + 1;
    end
    if(j ~= 48)
        Azcorr(rz) = corr(48, j);
        rz = rz + 1;
    end
end

figure();
subplot(3, 1, 1);
bar(linspace(1, 47, 47), Axcorr);
grid on;
xlabel('X index');
ylabel('\rho \Delta_R')
subplot(3, 1, 2);
bar(linspace(1, 47, 47), Aycorr);
grid on;
xlabel('X index');
ylabel('\rho \Delta_T')
subplot(3, 1, 3);
bar(linspace(1, 47, 47), Azcorr);
grid on;
xlabel('X index');
ylabel('\rho \Delta_N');
sgtitle('Coefficient correlation for position deviation');

% number of subplots
[p,~]=numSubplots(maxCount);
Nrow = p(1);
Ncol = p(2);

% plot error per iteration
figure();
lw = 1.5;
for c = 1:maxCount
    subplot(Nrow, Ncol, c);
    up = 3 * c;
    down = up - 2;
    errk = errj(down:up, :);
    plot(t, errk(1, :), LineWidth= lw);
    hold all;
    plot(t, errk(2, :), LineWidth= lw);
    plot(t, errk(3, :), LineWidth= lw);
    xlabel('TIME [s]');
    ylabel('position error [m]');
    grid on;
    title('Iteration ' + string(c));
    legend('R', 'T', 'N');
end
sgtitle('absolute position error per iteration');

% plot RMS coeff error
figure();
lw = 1.5;
nval = linspace(2, n_max, n_max-1);
for c = 1:maxCount
    subplot(Nrow, Ncol, c);
    if(maxCount ==  4) 
        hold on;
        semilogy(nval, coefErr_truth_noise_E12(c, 2:end), 'o --', LineWidth= lw);
    end
    semilogy(nval, CoefErr(c, 2:end), 'sq --', LineWidth= lw);
    xlabel('Harmonic degree');
    ylabel('error RMS');
    grid on;
    title('Iteration ' + string(c));
end
sgtitle('RMS coefficient error per iteration');
if(maxCount == 4)
    legend('No deviation', 'State deviation included');
else
    legend('ODEST');
end

% create table
tt = [nval', CoefErr(end, 2:end)', coefErr_truth_noise_E12(end, 2:end)'];
T = array2table(tt,...
    'VariableNames',{'Degree','SDIM','No deviation orbit'});
disp(T);



%% FUNCTIONS
function [CoefErr] = computeRMS_coefErr(n_max, Nc, Ns, X, Cnm, Snm)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = (X(1) - 5.2)^2;
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
    for n = 2:n_max
        CoefErr(n) = sqrt(CoefErr(n) / (2*n + 1));
    end
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

