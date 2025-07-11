clear;
clc;
close all;
format long g;

addpath('functions/');
addpath('../QGG_gravEstim/src/');
addpath('../QGG_gravEstim/data_files/');
addpath('../matlab_codes/GOCE_products/GOCE_L2b_MatlabReaders/data/');

set(0,'defaultAxesFontSize',16);

%%      EVALUATE NOISE PROJECTION
% Description: Undertand and evaluate noise projection in Null space
% method.
% Author: Sergio Coll
% Date: 05/21/25

% select error type
Planet   = "Earth";       % options: Earth / Bennu / Eros
error    = "attitude";    % options: position / attitude
saveData = 0;             % options: 0 / 1

sensitivity_analysis = 0;   % options: 0 /1
residuals_projection = 1;   % options: 0 /1
information_loss     = 0;   % options: 0 /1

[planetParams, poleParams, Kaula, r, Xtrue] = loadPlanet(Planet);
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 

W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
DEC = poleParams(4);

% Time options
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 1;
f   = 1/60;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);

% Gravity field parameters
[Cnm, Snm] = list2mat(n_max, Nc, Ns, Xtrue);

if(saveData == 0)
    r0 = [r;0;0];           % [ACI]
    v0 = [0;0;sqrt(GM/r)];  % [ACI]
    
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    PHI0 = reshape(eye(6,6), [36, 1]);
    [~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 4, GM, Re, normalized, ...
        W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
    ACI_ECEF = 0;
    rn = state_t(:, 1:3)';
    vn = state_t(:, 4:6)';
else
    load('Nov_L2position.mat');   % ECEF coordinates
    load('Nov_L2velocity.mat');   % ECEF coordinates
    load('Nov_L2ECEF2ITRF.mat');  % rotation matrix ECEF 2 ITRF
    load('Nov_L2ITRF2GRF.mat');   % rotation matrix ITRF 2 GRF
    
    [~, t2]  = quaternion2CDM(outPut);
    [~, t3]  = quaternion2CDM(outPut2);
    
    % Check time-points to make sure that all files are at the
    % same time.
    t1 = positions(:, 1);
    commonTimes = intersect(intersect(t1, t2), t3);
    
    [~, idx1] = ismember(commonTimes, t1);
    [~, idx2] = ismember(commonTimes, t2);
    [~, idx3] = ismember(commonTimes, t3);

    t = t1(idx1)';
    Nt = length(t);
    
% %     % WARNING: reducing data set
% %     t = t(1:round(Nt/10));
% %     Nt = length(t);

    rn_ECEF = positions(idx1, 2:end)';
    vn_ECEF = velocity(idx1, 2:end)';
    ACI_ECEF = quaternion2CDM(outPut(idx2, :));
    GRF_ACI  = quaternion2CDM(outPut2(idx3, :));

    state_ACI = rotate2ECI([rn_ECEF; vn_ECEF], ACI_ECEF, t);
    state_t(:, 1:3) = state_ACI(1:3, :)';
    state_t(:, 4:6) = state_ACI(4:6, :)';
end

% compute planet orientation parameters
[ACAF_ACI] = compute_planetAtt(poleParams, t, saveData, ACI_ECEF);

% compute noise sensitivity
if(sensitivity_analysis)
    [Hmat] = compute_Noise_Sensitivity(t, rn, vn, ACAF_ACI, Cnm,...
        Snm, Re, GM, n_max, normalized, error);
    for j = 1:1
        figure()
        y = smoothdata(Hmat(:, :, j), 1, 'movmean', 1);
        z = mean(Hmat(:, :, j));
        plot(t, y, 'LineWidth', 2)
        if(error == "position")
            legend('', '$S_{\sigma^2_{RT}} , \bar{S} = $' +string(z(2)), '$S_{\sigma^2_{RN}} , \bar{S} = $' + string(z(3)),...
                '$S_{\sigma^2_{TT}} , \bar{S} = $' + string(z(4)), '$S_{\sigma^2_{TN}} , \bar{S} = $' + string(z(5)), '$S_{\sigma^2_{NN}} , \bar{S} = $' + string(z(6)), 'Interpreter', 'latex');
        else
            legend('$S_{\sigma^2_{RR}} , \bar{S} = $' + string(z(1)), '$S_{\sigma^2_{RT}} , \bar{S} = $' +string(z(2)), '$S_{\sigma^2_{RN}} , \bar{S} = $' + string(z(3)),...
                '$S_{\sigma^2_{TT}} , \bar{S} = $' + string(z(4)), '$S_{\sigma^2_{TN}} , \bar{S} = $' + string(z(5)), '$S_{\sigma^2_{NN}} , \bar{S} = $' + string(z(6)), 'Interpreter', 'latex');
        end
    end
end

%  compute residual projecction
if(residuals_projection)
    [coeffs_vals, terms] = compute_coeff_proj(t, rn, vn, ACAF_ACI, ...
        Cnm, Snm, Re, GM, n_max, normalized, error);
    
    % plotting
    for j = 1:length(coeffs_vals(1, 1, :))
        figure()
        plot(t, coeffs_vals(:, :, j)', 'LineWidth', 2)
        if(length(terms) == 6)
            legend('dY_{RR}','dY_{RT}', 'dY_{RN}', 'dY_{TT}', 'dY_{TN}', 'dY_{NN}')
        else
            legend('dY_{RT}', 'dY_{RN}', 'dY_{TT}', 'dY_{TN}', 'dY_{NN}')
        end
        
        title('$dY \cdot \vec{v}_{' + string(j) + '}$', 'Interpreter','latex');
    end
end

% compute information loss after the projection
if(information_loss)
    [I_proj, I_org, I_proj_ideal] = compute_Inf_loss(t, rn, vn, ACAF_ACI,...
        Cnm, Snm, Re, GM, n_max, normalized, error);
    if(n_max <=2)
        plot_information_lowGravField(I_proj, I_org, I_proj_ideal, ...
            n_max, Nc, Ns, Ncs);
    else
        plot_information_highGravField(I_proj, I_org,...
            n_max, Nc, Ns)
    end
end



%% FUNCTIONS

function [] = plot_information_highGravField(I_proj, I_org, n_max, Nc, Ns)
    ratio = (1 - diag(I_proj./I_org)) * 100;
    [A, B] = list2mat(n_max, Nc, Ns, ratio);
    A(1,1) = ratio(1);

    J = ones(n_max+1, (n_max+1)*2) * NaN; S = fliplr(B);
    J(:, 1:n_max+1) = S;
    J(:, n_max+2:end) = A;

    figure;
    cmap = turbo(256);  % Original colormap with 256 colors
    im = imagesc(J);
    clim([1 100]);  % Force color scaling between 1% and 100%
    colormap(cmap);
    set(im, 'AlphaData', ~isnan(J)); % Make NaNs transparent
    hold on;

    c = colorbar;
    c.Ticks = [1, 20, 40, 60, 80, 100]; 
    c.TickLabels = {'1', '20', '40', '60', '80', '100 %'};
   
    % Set axis properties
    yticks(linspace(1, n_max + 1, 7));
    yticklabels(compose('%i', linspace(0, n_max, 7)));
    xticks(linspace(1, 2*n_max + 2, 7));
    xticklabels(compose('%i', linspace(-n_max, n_max, 7)));
    
    xlabel('Order');
    ylabel('Degree');
    hold off;
    title('Information loss due to projection');
end

function [] = plot_information_lowGravField(I_proj, I_org, I_proj_ideal, n_max, Nc, Ns, Ncs)
    % Plot results
    [num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
    
    % Account for 00 term
    num_C = num_C + 1;
    num_C = [1, num_C];
    str_C = ["C_{00}", str_C];
    
    figure()
    B = 1 - diag(I_proj./I_org); % information ratio
    subplot(1, 2, 1)
    plot(1:Nc, B(1:Nc), 'LineWidth', 2)
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    title('C_{nm}')
    subplot(1, 2, 2)
    plot(1:Ns, B(Nc+1:Ncs), 'LineWidth', 2)
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    title('S_{nm}')
    sgtitle('Information loss due to projection')
    
    figure()
    B = diag(I_proj./I_proj_ideal);           % projected information ratio
    subplot(1, 2, 1)
    semilogy(1:Nc, B(1:Nc), 'LineWidth', 2)
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    title('C_{nm}')
    subplot(1, 2, 2)
    semilogy(1:Ns, B(Nc+1:Ncs), 'LineWidth', 2)
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    title('S_{nm}')
    sgtitle('Projected information ratio w.r.t spherical case')
end

function [coeffs_vals, terms] = compute_coeff_proj(t, rn, vn, Rot, Cnm, Snm, Re, GM, n_max, normalized, error)
    syms y11 y12 y13 y22 y23 y33 real
    coeffs_vals = ones(length(t), 6, 3) * NaN;
    for k = 1:length(t)
        disp('Running .. ' + string(k/length(t)));
        
        % position and velcoity
        r  = rn(:, k); v = vn(:, k);
    
        % rotation to RTN frame
        ACI_RTN = RTN2ECI(r, v);
    
        % ACAF to ACI rotation matrix
        maxPos = 3 * k; minPos = maxPos- 2;
        ACAF_ACI = Rot(minPos:maxPos, :);
    
        % rotation to body frame
        ACAF_BODY = ACAF_ACI * ACI_RTN;
% %         ACAF_RTN = ACAF_ACI;
    
        % position / attitude partials & null space
        if(error == "position")
            [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
                 ACAF_ACI, ACAF_BODY);
             C = null([Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
             dY = [y12 y13 y22 y23 y33]';

% %              C = null([Hp(1,:);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
% %              dY = [y11 y12 y13 y22 y23 y33]';
        elseif(error == "attitude")
            [Hp] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_BODY);
            C = null([Hp(1,:); Hp(2,:);Hp(3,:);Hp(5, :); Hp(6, :); Hp(9, :)]');
            dY = [y11 y12 y13 y22 y23 y33]';
        else
            warning('Select a correct mode for the error type')
            break;
        end
        Nterms  = length(dY);
        Nspace  = Nterms - 3;

        % meas residual proyection
        dY_proj  = C' * dY;
    
        % get the coefficients
        for j = 1:Nspace
            [coeffs_vals(k, 1:Nterms, j), terms] = coeffs(dY_proj(j), dY);
        end
    end
end

function [Hmat] = compute_Noise_Sensitivity(t, rn, vn, Rot, Cnm, Snm, Re, GM, n_max, normalized, error)
    % noise values
    Hmat = ones(length(t), 6 , 3) * NaN;
    
    % create noise values
    syms v11 v12 v13 v21 v22 v23 v31 v32 v33 positive
    R = diag([v11,v12,v13,v22,v23,v33]);

    for k = 1:length(t)
        disp('Running .. ' + string(k/length(t)));
        
        % position and velcoity
        r  = rn(:, k); v = vn(:, k);
    
        % rotation to RTN frame
        ACI_RTN = RTN2ECI(r, v);
    
        % ACAF to ACI rotation matrix
        maxPos = 3 * k; minPos = maxPos- 2;
        ACAF_ACI = Rot(minPos:maxPos, :);
    
        % rotation to body frame
        ACAF_RTN = ACAF_ACI * ACI_RTN;
% %         ACAF_RTN = ACAF_ACI;
        if(error == "position")
            [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
            ACAF_ACI, ACAF_RTN); 
            C = null([Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
            R = diag([v12,v13,v22,v23,v33]);
        elseif(error == "attitude")
            [Hp] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_RTN);
            C = null([Hp(1, :); Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
            R = diag([v11, v12,v13,v22,v23,v33]);
        else
            warning('Select a correct error mode...');
        end
        % noise projection
        r = C' * R * C;
        d = det(r);

        for i = 1:length(d)
            % compute the sensitivity matrix
            J11 = diff(d(i), v11);
            J12 = diff(d(i), v12);
            J13 = diff(d(i), v13);
            J22 = diff(d(i), v22);
            J23 = diff(d(i), v23);
            J33 = diff(d(i), v33);
        
            % evaluate around a nominal sensitivity value
            val = 1;
            H11 = subs(J11, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H12 = subs(J12, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H13 = subs(J13, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H22 = subs(J22, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H23 = subs(J23, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H33 = subs(J33, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
        
            Hmat(k, :, i)  = [H11, H12, H13, H22, H23, H33];
        end
    end
end

function [I_proj, I_org, I_proj_ideal] = compute_Inf_loss(t, rn, vn, Rot, Cnm, Snm, Re, GM, n_max, normalized, error)
    I_proj = 0; I_org = 0; I_proj_ideal = 0;
    sigma_RR = 1; sigma_RT = 1; sigma_RN = 1;
    sigma_TT = 1; sigma_TN = 1; sigma_NN = 1;
    for k = 1:length(t)
        disp('Running .. ' + string(k/length(t)));
        
        % position and velcoity
        r  = rn(:, k); v = vn(:, k);
    
        % rotation to RTN frame
        ACI_RTN = RTN2ECI(r, v);
    
        % ACAF to ACI rotation matrix
        maxPos = 3 * k; minPos = maxPos- 2;
        ACAF_ACI = Rot(minPos:maxPos, :);
    
        % rotation to body frame
        ACAF_BODY = ACAF_ACI * ACI_RTN;
        ACI_BODY  = ACI_RTN;
% %         ACAF_BODY = ACAF_ACI;
% %         ACI_BODY = eye(3,3);

        % compute gravity field partials
        planetParams = [GM, Re, n_max, normalized];
        [~, ~, H] = gradiometer_meas(t(k) ,planetParams, ACAF_ACI, [r', v'], ...
                zeros(9, length(t)), Cnm, Snm, ACI_BODY);
    
        % position partials & null space
        if(error == "position")
            [Hp] = compute_posPartials(0, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
                 ACAF_ACI, ACAF_BODY);
             C = null([Hp(2,:);Hp(3,:);Hp(5,:);Hp(6, :);Hp(9, :)]');
             R = diag([sigma_RT, sigma_RN, sigma_TT,sigma_TN,sigma_NN].^2);
             R_id = eye(5,5);
             h = [H(4,:);H(7,:);H(5,:);H(8, :);H(9,:)];
             
% %              C  = null([Hp(3,:);Hp(5,:);Hp(6,:);Hp(9, :)]');
% %              R = diag([sigma_RN,sigma_TT, sigma_TN,sigma_NN].^2);
% %              R_id = eye(4,4);
% %              h = [H(7,:);H(5,:);H(8,:);H(9, :)];

             h_org = [H(1, :);H(4,:);H(7,:);H(5, :);H(8, :);H(9, :)];
             R_org = diag([sigma_RR,sigma_RT, sigma_RN, sigma_TT, sigma_TN, sigma_NN].^2);
% %              h_org = h;
% %              R_org = R;
             
        elseif(error == "attitude")
            [Hp] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_BODY);
            C = null([Hp(1,:); Hp(2,:);Hp(3,:);Hp(5, :); Hp(6, :); Hp(9, :)]');
            R = diag([sigma_RR, sigma_RT, sigma_RN, sigma_TT, sigma_TN, sigma_NN].^2);
            R_id = eye(6,6);
            h = [H(1, :);H(4,:);H(7,:);H(5, :);H(8, :); H(9, :)];


% %             C = null([Hp(2,:);Hp(3,:);Hp(5, :);Hp(6,:);Hp(9,:)]');
% %             R = diag([sigma_RT,sigma_RN,sigma_TT, sigma_TN, sigma_NN].^2);
% %             R_id = eye(5,5);
% %             h = [H(4,:);H(7,:);H(5,:);H(8,:);H(9, :)];
% % 
            h_org = [H(1, :);H(4,:);H(7,:);H(5, :);H(8, :); H(9, :)];
            R_org = diag([sigma_RR,sigma_RT, sigma_RN, sigma_TT, sigma_TN, sigma_NN].^2);
% %             h_org = h;
% %             R_org = R;
        else
            warning('Select a correct mode for the error type')
            break;
        end

        % project matrices
        h_proj = C' * h;
        r = C' * R * C;
        r_ideal = C' * R_id * C;

        % compute information
        I_org = I_org + h_org' * inv(R_org) * h_org;
        I_proj = I_proj + h_proj' * inv(r) * h_proj;
        I_proj_ideal = I_proj_ideal + h_proj' * inv(r_ideal) * h_proj;
    end
end