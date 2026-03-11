clear; clc;
close all;

addpath("../../att_errors_res_NSM/functions/");
addpath("../../regress_reg/");
addpath("../../GG_apriori/");
%% MODIFY REGRESS FILE

% input regress files
file_inReg_XX  = "goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";
file_inReg_XY  = "goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";
file_inReg_XZ  = "goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";
file_inReg_YY  = "goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";
file_inReg_YZ  = "goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";
file_inReg_ZZ  = "goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";

file_outReg    = "test_goce_YY_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg";

% input observation file
file_inObs     = "goce_eggreg_2008-08-01_RL5061_RL05surpv6.980.001.ggr";
GG             = read_GG_obs(file_inObs);

% Read regress partials + residuals + sigma
H_XX           = read_reg_files_partials(file_inReg_XX);
H_XY           = read_reg_files_partials(file_inReg_XY);
H_XZ           = read_reg_files_partials(file_inReg_XZ);
H_YY           = read_reg_files_partials(file_inReg_YY);
H_YZ           = read_reg_files_partials(file_inReg_YZ);
H_ZZ           = read_reg_files_partials(file_inReg_ZZ);

% projected partials
H_XX_proj = H_XX;
H_XY_proj = H_XY;
H_XZ_proj = H_XZ;
H_YY_proj = H_YY;
H_YZ_proj = H_YZ;
H_ZZ_proj = H_ZZ;

Nobs = length(H_XX(:, 1)); % Number of observations
mask = [1;1;1;0;1;1;0;0;1];
for k = 1:Nobs
    % disp( string(k/Nobs*100));
    % partials
    h = [H_XX(k, 1:end-2); H_XY(k, 1:end-2); H_XZ(k, 1:end-2);...
         H_YY(k, 1:end-2); H_YZ(k, 1:end-2); H_ZZ(k, 1:end-2)];
    
    % residuals
    rs= [H_XX(k, end-1); H_XY(k, end-1); H_XZ(k, end-1);...
         H_YY(k, end-1); H_YZ(k, end-1); H_ZZ(k, end-1)];

    % observations
    g   = GG(k, 2:end)';
    obs = [g(1);g(2);g(3);g(2);g(4);g(5);g(3);g(5);g(6)];
    
    % orientation error partials
    [Hrot] = compute_rotPartials_analy(obs, eye(3));
    hrot   = Hrot(logical(mask),:);
    [s, v, d]  = svd(hrot');
    V_rot  = d(:, 3:6);
    P_rot  = V_rot * V_rot';

    % project partials and residuals
    h_proj  = P_rot * h;
    rs_proj = P_rot * rs;

    % write projected partials
    H_XX_proj(k, 1:end-1) = [h_proj(1, :), rs_proj(1, :)];
    H_XY_proj(k, 1:end-1) = [h_proj(2, :), rs_proj(2, :)];
    H_XZ_proj(k, 1:end-1) = [h_proj(3, :), rs_proj(3, :)];
    H_YY_proj(k, 1:end-1) = [h_proj(4, :), rs_proj(4, :)];
    H_YZ_proj(k, 1:end-1) = [h_proj(5, :), rs_proj(5, :)];
    H_ZZ_proj(k, 1:end-1) = [h_proj(6, :), rs_proj(6, :)];
end



parCols  = 1:length(H(1,:))-2;
newVals  = H.*30;


info = reg_update_partials_allRows_v2(file_inReg,file_outReg,parCols,newVals);
disp(info)

H_new    = read_reg_files_partials(file_outReg);