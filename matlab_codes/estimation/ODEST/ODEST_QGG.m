clc;
close all;
clear;
format long g;

%%                          ODEST FUNCTION TEST
% Description: test script for the ODEST algorithm using data processed 
% with the QGG code.

% INPUT
addpath('../functions/');
addpath('../dataFiles/');
addpath('../../../QGG_gravEstim/data_files/') ;

dOrbit = readtable("orbitData.txt");
dAcc = readtable('accData.txt');
dataOrbit = table2array(readtable("orbitData.txt"));
dataMeas = table2array(readtable('accData.txt'));

set(0,'defaultAxesFontSize',16);

% Simualtion parameters
n_max = 6;
normalized = 0;
option = "required";       % noise option: required / obtained
sigma = 6.3E-13;           % Noise level

% Harmonics values. Truth (matrix form)
path = "HARMCOEFS_BENNU_OSIRIS_0.txt";
[Cnm_t, Snm_t, Re] = readCoeff(path);

% Bennu parameters
rho = 1250;
G_constant = 6.6743E-11;
% GM = rho * 4/3 * pi * (Re)^3 * G_constant;
GM = 5.2;
n = sqrt(GM / 1000^3);
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% Harmonic values. A priori (vector form)
[Nc, Ns, Ncs] = count_num_coeff(n_max);
[Cnm_ap, Snm_ap] = getkaula(n_max, Nc, Ns, normalized);

% number of points
t = dataOrbit(:, 1);
At = t(2) - t(1);
Nt_max = length(t);

% truth position
rt_ACI = dataOrbit(:, 2:4)';
vt_ACI = [dataOrbit(:, 8), dataOrbit(:, 9), dataOrbit(:, 10)]'; 

% Nominal orbit. ODE 113 integrator
X0 = [rt_ACI(:, 1); vt_ACI(:, 1)];
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, Cnm_t, Snm_t, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, X0, options);
state = state';

% print nominal vs truth
figure();
for k = 1:3
    subplot(3, 1, k)
    plot(t./86400, state(k, :), t./86400, rt_ACI(k, :), 'LineWidth', 1.5);
    grid on;
    legend('Nominal', 'Truth');
    xlabel('TIME [days]')
    if(k == 1)
        ylabel('I [m]')
    elseif(k == 2)
        ylabel('J [m]')
    else
        ylabel('K [m]')
    end
end
sgtitle('Truth vs nominal. Inertial frame, ACI');

figure();
for k = 1:3
    err = state(k, :) - rt_ACI(k, :);
    subplot(3, 1, k)
    plot(t./86400, err, 'LineWidth', 3, ...
        'Color', 'red');
    grid on;
    xlabel('TIME [days]')
    if(k == 1)
        ylabel('I [m]')
    elseif(k == 2)
        ylabel('J [m]')
    else
        ylabel('K [m]')
    end
end
sgtitle('Absolute error (Truth - nominal). Inertial frame, ACI');

% nominal orbit. ACI
rr_ACI = [state(1, :); state(2,:); state(3, :)];
vrr_ACI = [state(4, :); state(5,:); state(6, :)];
rr_ACAF = zeros(3, Nt_max);
rr_RTN = zeros(3, Nt_max);
rt_RTN = zeros(3, Nt_max);

% rotation matrix
ACAF_RTN_mat = zeros(3*Nt_max, 3);
RTN_ACI_mat = zeros(3*Nt_max, 3);
ACAF_ACI_mat = zeros(3*Nt_max, 3);

% bias values correction
[B, D] = getBias(option);

% measurements vector. Truth
Y_o = zeros(9, Nt_max);
Y_oc = zeros(9, Nt_max);
for k = 1:Nt_max
    [ACI_RTN] = RTN2ECI(rt_ACI(:, k), vt_ACI(:, k));
    RTN_ACI = ACI_RTN';
    up = 3 * k;
    down = up - 2;
    
    rt_RTN(:, k) = RTN_ACI * rt_ACI(:, k);

    Wt = W0 + W * t(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    ACAF_ACI_mat(down:up, :) = ACAF_ACI;

    y = -[dataMeas(k, 2), dataMeas(k, 3), dataMeas(k, 4);...
    dataMeas(k, 5), dataMeas(k, 6), dataMeas(k, 7);...
    dataMeas(k, 8), dataMeas(k, 9), dataMeas(k, 10)] + B + ...
    (k - 1) * At * D; % WARNING. I had to change the sign

    Y_o(:, k) = reshape(y, [9 ,1]);
end

% compute rotation matrices and position at nominal
for k =1:Nt_max
    up = 3*k;
    down = up - 2;

    [ACI_RTN] = RTN2ECI(rr_ACI(:, k), vrr_ACI(:, k));
    RTN_ACI = ACI_RTN';

    ACAF_RTN = ACAF_ACI_mat(down:up, :) * ACI_RTN;

    % save rotation matrices
    ACAF_RTN_mat(down:up, :) = ACAF_RTN;
    RTN_ACI_mat(down:up, :) = RTN_ACI;

    % position at nominal
    rr_ACAF(:, k) = ACAF_ACI_mat(down:up, :) * rr_ACI(:, k);
    rr_RTN(:, k) = RTN_ACI * rr_ACI(:, k);
end


% inital values (vector form)
X = [Cnm_ap; Snm_ap];
[Cnm_ap, Snm_ap] = list2mat(n_max, Nc, Ns, X);

% ODEST 
count = 1;
maxCount = 4;
R = [2*sigma^2, 0;...
     0, sigma^2];

Aj = zeros(3*maxCount, Nt_max);
errj = zeros(3*maxCount, Nt_max);
CoefErr = zeros(maxCount, n_max);
Xx = 0;
Nx = Ncs;
while(count <= maxCount)
    disp('      Iteration: ' + string(count));
    % Batch values
    P0 = eye(Nx - 1) * 1E1;
    if(count == 1)
        X0 = zeros(Nx - 1, 1);
    else
        X0 = Xx;
    end
    Axp = inv(P0);
    Nxp = P0\-X0;
    
    % time loop
    for k = 1:Nt_max
        % read rotation matrices
        up = 3*k;
        down = up - 2;
        ACAF_RTN = ACAF_RTN_mat(down:up, :);
        
        % RTN ACI
        RTN_ACI = RTN_ACI_mat(down:up, :);

        % ACAF ACI
        ACAF_ACI = ACAF_ACI_mat(down:up, :);

        % compute position
        rr_RTN(:,  k) = RTN_ACI * rr_ACI(:, k);

        % compute perturbed meas. RTN frame
        [~, ~, ddU] = potentialGradient_nm(Cnm_ap, Snm_ap, n_max, ...
                                               rr_ACAF(:, k), Re, GM, ...
                                               normalized);
        
        Yc = ACAF_RTN'* ddU * ACAF_RTN;
        Yc = reshape(Yc, [9, 1]);

        % rotate meas to RTN. Nominal
        y = reshape(Y_o(:, k), [3,3]);
        y = reshape(RTN_ACI * y * RTN_ACI', [9,1]);
    
        % H perturbed. ENU
        [Hp] = potentialGradient_Cnm(n_max, rr_ACAF(:, k), Re, GM,...
            ACAF_RTN', normalized);
    
        % new meas
        Y = y - Yc;

        % pre-weighting meas
        Y = [Y(5) - Y(9); Y(6)];
        Hp = [Hp(5, 2:end) - Hp(9, 2:end); Hp(6, 2:end)];
        
        Axp = Axp + (Hp' * inv(R) * Hp);
        Nxp = Nxp + (Hp' * inv(R) * Y);
    end
    
    % Harmonic correction
    XNOT = Axp\Nxp;
    
    % update state
    X(2:end) = X(2:end) + XNOT;
    [Cnm_ap, Snm_ap] = list2mat(n_max, Nc, Ns, X);

    % ODEST algorithm
    [AR, AN, AT] = ODEST(Nt_max, n_max, Re, GM,  ...
        rr_ACI, rr_ACAF, ACAF_RTN_mat, RTN_ACI_mat, Cnm_ap, ...
        Snm_ap, Y_o);

% %     [AR, AT, AN] = ODEST_arc(Nt_max, n_max, Y_o, Re, GM, ...
% %             rr_ACAF, Cnm_updt, Snm_updt, ACAF_RTN_mat, RTN_ACI_mat, sigma);

    A = [AR; AT; AN];
    up = count*3;
    down = up - 2;
    Aj(down:up, :) = A;
    errj(down:up, :) = rt_ACI - rr_ACI;

    % update position
    for k = 1:Nt_max
        up = 3*k;
        down = up -2;

        A_k = A(:, k);
        rr_ACAF(:, k) = rr_ACAF(:, k) - ACAF_RTN_mat(down:up, :) * A_k;
        rr_ACI(:, k) = ACAF_ACI_mat(down:up, :)' * rr_ACAF(:, k);
    end
    
    % update rotation matrices. Nominal
     vrr_ACI = computeVel(rr_ACI, t);
    for k =1:Nt_max
        up = 3*k;
        down = up - 2;
    
        [ACI_RTN] = RTN2ECI(rr_ACI(:, k), vrr_ACI(:, k));
        RTN_ACI = ACI_RTN';
    
        ACAF_RTN = ACAF_ACI_mat(down:up, :) * ACI_RTN;
    
        % save rotation matrices
        ACAF_RTN_mat(down:up, :) = ACAF_RTN;
        RTN_ACI_mat(down:up, :) = RTN_ACI;
    end
    
    Xx = Xx + XNOT;
    if(n_max == 6)
        CoefErr(count, :) = computeRMS_coefErr(n_max, 26, 20, X, Cnm_t, Snm_t);
    elseif(n_max == 3)
        CoefErr(count, :) = computeRMS_coefErr(n_max, 8, 5, X, Cnm_t, Snm_t);
    end

    % update counter
    count = count + 1;
end

%% PLOTS

% plot correction
[p,~]=numSubplots(maxCount);
Nrow = p(1);
Ncol = p(2);

figure();
lw = 1.5;
for c = 1:maxCount
    subplot(Nrow, Ncol, c);
    up = 3 * c;
    down = up - 2;
    Ak = Aj(down:up, :);
    plot(t, Ak(1, :), LineWidth= lw);
    hold all;
    plot(t, Ak(2, :), LineWidth= lw);
    plot(t, Ak(3, :), LineWidth= lw);
    xlabel('TIME [s]');
    ylabel('deviation prediction [m]');
    grid on;
    title('Iteration ' + string(c));
    if(c == 1)
        legend('\Delta_R', '\Delta_T', '\Delta_N');
    end
end
sgtitle('Postion deviation prediction per iteration');

% plot error per iteration
figure();
lw = 3;
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
    if(c == 1)
        legend('I', 'J', 'K');
    end
end
sgtitle('absolute position error per iteration. ACI frame');

% plot RMS coeff error
figure();
lw = 3;
nval = linspace(2, n_max, n_max-1);
for c = 1:maxCount
    subplot(Nrow, Ncol, c);
    semilogy(nval, CoefErr(c, 2:end), 'sq --', 'Color', 'b', LineWidth=lw);
    xlabel('Harmonic degree');
    ylabel('error RMS');
    grid on;
    title('Iteration ' + string(c));
end
sgtitle('RMS coefficient error per iteration');
if(maxCount == 4)
    legend('ODEST');
end

% absolute error deviation predicton
figure();
err = rt_RTN - rr_RTN;
errR = (rt_RTN - rr_RTN)./rt_RTN;
subplot(3, 1, 1);
plot(t, err(1, :) , LineWidth=1.5);
grid on;
xlabel('TIME [s]')
ylabel('R error');

subplot(3, 1, 2);
plot(t, err(2, :),  LineWidth=1.5);
grid on;
xlabel('TIME [s]')
ylabel('T error');

subplot(3, 1, 3);
plot(t, err(3, :), LineWidth=1.5);
grid on;
xlabel('TIME [s]')
ylabel('N error');
sgtitle('position absolute error using ODEST [m]')

% create table
tt = [nval', CoefErr(end, 2:end)'];
T = array2table(tt,...
    'VariableNames',{'Degree','ODEST: ON'});
disp(T);


%% FUNCTION

function [AR, AN, AT] = ODEST(Nt_max, n_max, Re, GM,  ...
    rr_ACI, rr_ACAF, ACAF_RTN_mat, RTN_ACI_mat, Cnm_updt, Snm_updt, Y_o)

    % output values. Position deviation
    AR = zeros(1, Nt_max);
    AT = zeros(1, Nt_max);
    AN = zeros(1, Nt_max);

    for k = 1:Nt_max
        % nominal radius
        rp = vecnorm(rr_ACI(:, k));
    
        % perturbed potential
        up = 3 * k;
        down = up - 2;
        R = ACAF_RTN_mat(down:up, :);
        RTN_ACI = RTN_ACI_mat(down:up, :);
        
        % using corrected harmonics
        [~, ~, delta_ddU] = potentialGradient_nm(Cnm_updt, Snm_updt, n_max, ...
                                            rr_ACAF(:, k), Re, GM, 0);
        delta_ddU = R' * delta_ddU * R;
        Y = reshape(delta_ddU, [9, 1]);
        
        y = reshape(Y_o(:, k), [3,3]);
        y = reshape(RTN_ACI * y * RTN_ACI', [9, 1]);
        Y_oc_new = y - Y;
        
        % position deviation
        AR(k) = rp^4 / (6*GM) *  (Y_oc_new(1));
        AN(k) = -rp^4 / (3*GM) * (Y_oc_new(3));
        AT(k) = -rp^4 / (3*GM) * (Y_oc_new(4));
    end

    % signal smoothing
% %     polinomial_order = 3;
% %     window_points = 2001;
% %     AR = sgolayfilt(AR, polinomial_order, window_points);
% %     AT = sgolayfilt(AT, polinomial_order, window_points);
% %     AN = sgolayfilt(AN, polinomial_order, window_points);
end

function [AR, AT, AN] = ODEST_arc(Nt_max, n_max, Y_o, Re, GM, ...
    rr_ACAF, Cnm_updt, Snm_updt, ACAF_RTN_mat, RTN_ACI_mat,sigma)
    % Description: ODEST algorithm using the arc method discretization.

    Np_arc = 5;           % Number points per arc

    % correction matrices
    AR = zeros(1, Nt_max);
    AT = zeros(1, Nt_max);
    AN = zeros(1, Nt_max);

    % compute correction in the arc
    N0 = 1;
    Nf = Np_arc;
    while N0 <= Nt_max
        Ax = 0;
        Nx = 0;
        W = sigma^2;

        for k = N0:Nf
            % radius vec
            rp = vecnorm(rr_ACAF(:, k));

            % rotation matrix
            up = 3 * k;
            down = up - 2;
            ACAF_RTN = ACAF_RTN_mat(down:up, :);
            RTN_ACI = RTN_ACI_mat(down:up, :);
            
            % compute value
            [~, ~, delta_ddU] = potentialGradient_nm(Cnm_updt, Snm_updt, n_max, ...
                                                rr_ACAF(:, k), Re, GM, 0);
            delta_ddU = ACAF_RTN' * delta_ddU * ACAF_RTN;
            Yc = reshape(delta_ddU, [9, 1]);

            yo = reshape(Y_o(:, k), [3,3]);
            yo = reshape(RTN_ACI * yo * RTN_ACI', [9, 1]);
            
            Y_oc_AR = yo(1) - Yc(1);
            Y_oc_AT = yo(3) - Yc(3);
            Y_oc_AN = yo(4) - Yc(4);

            H_AR = - 6 * GM / (rp^4);
            H_AT = 3 * GM / (rp^4);
            H_AN = H_AT;

            deltaY = [Y_oc_AR; Y_oc_AT; Y_oc_AN];
            Hi = [H_AR, 0, 0;...
                0, H_AT, 0;...
                0, 0, H_AN];

            Ri_inv = inv(W);

            Ax = Ax + (Hi' * Ri_inv * Hi);
            Nx = (Hi' * Ri_inv * deltaY) + Nx;
        end
        % solve LLS
        A  = -Ax\Nx;

        % compute corrections. Constant in the arc
        AR(N0:Nf) = A(1);
        AN(N0:Nf) = A(2);
        AT(N0:Nf) = A(3);
    
        % update arc
        N0 = Nf + 1;
        Nf = N0 + Np_arc;
        if(Nf > Nt_max)
            Nf = Nt_max;
        end
    end
end

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

function [p,n] = numSubplots(n)
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

function [B, D] = getBias(option)
    if(option == "obtained")
        % Obtained
        B = [-1.41000139962798e-07;...
            4.42692393515439e-09;...
            8.21383297049534e-09;...
            -3.18001688021536e-07;...
            -9.74999439415953e-06;...
            -5.65000358162869e-08];
        
        D = [ 3.78130353377026e-17;...
            8.11853407415446e-17;...
            -6.23503707619306e-18;...
            1.13081945673253e-15;...
            5.41557860429089e-15;...
            1.96147289639821e-17];
    elseif(option == "required")
        % Required
        B = [-282,8.88,16.42;1848,-636,-19500;5,9699,-113] ...
        * 1E-9/2;  % [1/s^2]

        % drift value. Truth
        D = [0.0065,0.0103, 0.0001; 3.123,0.195,0.9374;0.001, 2.98,0.0033] ...
            * 1E-9/(2*86400); % [1/s^2/s]
    end
end
