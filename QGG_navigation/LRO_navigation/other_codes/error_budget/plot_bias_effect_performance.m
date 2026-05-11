clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);

%% PLOT INLFUENCE OF BIAS IN POSITION PERFORMANCE
% Author: Sergio Coll-Ibars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = "/Users/sergiocollibars/OneDrive - UCB-O365/Kepler_codes/navigation_LRO/results/";

% % % set 1
% % file   = ["output_orbit_LRO_fc_1e-5_1mE.mat", "output_orbit_LRO_fc_1e-4_1mE.mat", ...
% %     "output_orbit_LRO_fc_1e-3_1mE.mat", "output_orbit_LRO_fc_1e-2_1mE.mat"];
% % fileNB   = ["NB_MOON_orbit_LRO_fc_1e-5_1mE.mat", "NB_MOON_orbit_LRO_fc_1e-4_1mE.mat", ...
% %     "NB_MOON_orbit_LRO_fc_1e-3_1mE.mat", "NB_MOON_orbit_LRO_fc_1e-2_1mE.mat"];
% % fileT  = "time_orbit_LRO_fc_1e-5_1mE.mat";

% % % set 2
% % file   = ["output_orbit_LRO_fc_1e-5_100mE.mat", "output_orbit_LRO_fc_1e-4_100mE.mat", ...
% %     "output_orbit_LRO_fc_1e-3_100mE.mat", "output_orbit_LRO_fc_1e-2_100mE.mat"];
% % fileNB   = ["NB_MOON_orbit_LRO_fc_1e-5_100mE.mat", "NB_MOON_orbit_LRO_fc_1e-4_100mE.mat", ...
% %     "NB_MOON_orbit_LRO_fc_1e-3_100mE.mat", "NB_MOON_orbit_LRO_fc_1e-2_100mE.mat"];
% % fileT  = "time_orbit_LRO_fc_1e-5_100mE.mat";


% set 3
file   = ["output_orbit_LRO_fc_1e-5_1e-2mE.mat", "output_orbit_LRO_fc_1e-4_1e-2mE.mat", ...
    "output_orbit_LRO_fc_1e-3_1e-2mE.mat", "output_orbit_LRO_fc_1e-2_1e-2mE.mat"];
fileNB   = ["NB_MOON_orbit_LRO_fc_1e-5_1e-2mE.mat", "NB_MOON_orbit_LRO_fc_1e-4_1e-2mE.mat", ...
    "NB_MOON_orbit_LRO_fc_1e-3_1e-2mE.mat", "NB_MOON_orbit_LRO_fc_1e-2_1e-2mE.mat"];
fileT  = "time_orbit_LRO_fc_1e-5_1e-2mE.mat";


load(folder + fileT);

T      = 2 * 3600;
epoch = time - time(1);
v = round(epoch./T);
maxorbs = max(v); prev = maxorbs - 1;
idx = find(v == maxorbs | v == prev);

fVec = [1E-5, 1E-4, 1E-3, 1E-2];
rmsVal = nan(3, length(fVec));
maxVal = nan(3, length(fVec));
minVal = nan(3, length(fVec));
for n = 1:length(file)
    load(folder+file(n));
    load(folder+fileNB(n));
   
    data   = mC_struct.sim1;
    Nt     = length(data.Xf(1,:));

    BN_mat = compute_orientation_SC(time, data.Xf(1:6,:)', "RTN", NB_MOON_mat);
    
    sigmaP = nan(3, Nt);
    for k = 1:Nt
        p = reshape(data.Pf(k, :), [18, 18]);
        maxInd = 3*k;
        minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);

        A   = blkdiag(BN, BN, eye(6),eye(6));
        P_B = A * p * A';
        std = sqrt(diag(P_B));
        sigmaP(:, k) = std(1:3);
    end
    rmsVal(:, n) = rms(sigmaP(:, idx)')';
    maxVal(:, n) = max(sigmaP(:, idx)')';
    minVal(:,n)  = min(sigmaP(:, idx)')';

% %     figure()
% %     plot(time(idx), sigmaP(:, idx));
end

figure(); tt = ["Radial", "Tangential", "Normal"];
clr = 'r';
for k = 1:3
    subplot(1, 3, k)
    loglog(fVec, rmsVal(k, :), 'Marker','x', 'LineStyle','-', 'Color', clr, 'LineWidth', 1.6);
    hold all;
    loglog(fVec, maxVal(k, :), 'Marker','o', 'LineStyle','none', 'Color', clr, 'MarkerFaceColor',clr);
    loglog(fVec, minVal(k, :), 'Marker','o', 'LineStyle','none', 'Color', clr, 'MarkerFaceColor',clr);
    grid on; xlabel('$f_c$', 'Interpreter','latex'); ylabel('[m]');
    ylim([1e-3 1e4]); title(tt(k));
end
sgtitle('Position 1\sigma uncertainty')