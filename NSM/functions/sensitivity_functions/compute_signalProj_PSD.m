function [] = compute_signalProj_PSD(planetParams, RotPlanet, B_ACI_mat, ...
        Xp, t, angVel, rn, vn, Iner, mask, Ytrue, att_err, attitude, angAcc)
     % planet variables
     n_max = planetParams(3);

     % variables number
    [Nc, Ns, ~] = count_num_coeff(n_max); 
     Nt          = length(t);

     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    % Projected signal 
    m = length(Ytrue(:, 1));
    r = 6;

    if(m > 6)
        Yproj = ones(m-r, Nt)*NaN;
    else
        Yproj = ones(1, Nt)*NaN;
    end
    Yproj2 = Ytrue.*0;
    res_att = Yproj;

    fprintf('       Progress:    0%%');  % Initial message
    for j = 1:Nt-2
        % Planet orientation
        maxPos = 3*j; minPos = maxPos - 2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);

        B_ACI = B_ACI_mat(minPos:maxPos, :);
   
        % Null space method correcting for attitude
        [Y_ACI, ~, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                zeros(9, Nt), Cp, Sp);
        [Yc] = add_angularComponents(Y_ACI, attitude(:, j), zeros(3, Nt), angVel(:, j).*0,...
                angAcc(:, j).*0);    

        % compute attitude partials. Nominal body frame
        [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
        [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
        Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];

        % apply maks to partials
        hrot = Hrot(logical(mask), :);
        
        % compute NS
        C = null(hrot');
% %         [~,~,d] = svd(hrot');
% %         C = d(:, 6);

        % projected signal
% %         Yproj(:, j) = (C')* Ytrue(:, j);
        Yproj(:, j) = (C')* Yc(logical(mask));
% %         Yproj2(:, j) = (C*C')* Ytrue(:, j);
        Yproj2(:, j) = (C*C')* Yc(logical(mask));

        % Attitude error residual signal
        res_att(:, j) = C' * hrot * att_err(:, j);

% %         Ytrue(:, j) = Yc(logical(mask));

         % Update every ~5% (optional)
        if mod(j, round(Nt/20)) == 0 || j == Nt-2
            fprintf('\b\b\b\b%3d%%', round(100 * j / Nt));
        end
    end

    % compute the PSD
    fs = 1/10;                 % 1 sample every 10 seconds

    window = hanning(512);
    nWindow = length(window);
    noverlap = 256;  % 75%
    nfft = 512*2;

    % Initialize
    [~, f] = pwelch(Ytrue(1,1:end-2), window, noverlap, nfft, fs);  % Frequency vector
    
    P_grav_total = zeros(nWindow+1, 1); P_grav_stored = zeros(nWindow+1, m);
    P_nsm_total  = zeros(nWindow+1, 1); P_nsm_stored = zeros(nWindow+1, r);
    P_res_stored = P_nsm_stored;

    TF = zeros(nWindow+1, m);
    for i = 1:m
        % Gravity-only signal
        [P_grav, ~] = pwelch(Ytrue(i,1:end-2), window, noverlap, nfft, fs);
        [P_nsm2, ~] = pwelch(Yproj2(i,1:end-2), window, noverlap, nfft, fs);

        P_grav_total = P_grav_total + P_grav;
        P_grav_stored(:, i) = P_grav;
        TF(:, i) = P_nsm2./P_grav;
% %         TF(:, i) = P_nsm2;
    end
    P_grav_total = P_grav_total * (1/m);

    for i = 1:length(Yproj(:, 1))
        % NSM-projected gravity
        [P_nsm, ~] = pwelch(Yproj(i,1:end-2), window, noverlap, nfft, fs);
        [P_res, ~] = pwelch(res_att(i,1:end-2), window, noverlap, nfft, fs);
        P_nsm_total = P_nsm_total + P_nsm;
        P_nsm_stored(:, i) = P_nsm;
        P_res_stored(:, i) = P_res;
    end
    P_nsm_total = P_nsm_total * (1/r);

    E_ratio = P_nsm_stored./ P_res_stored;

    figure();
    semilogx(f, 10*log10(TF));
    hold on;
    xlabel('Frequency [Hz]');
    ylabel('NSM Transfer Function [dB]');
    title('NSM Gravity Signal Transfer Function');
    grid on;

    figure(); % shows how much energy is contained in the original signal
    for j = 1:m
        semilogx(f, 10*log10(P_grav_stored(:, j)));
        hold on;
    end
    xlabel('Frequency [Hz]');
    ylabel('[dB]');
    title('Original Gravity Signal PSD');
    grid on;

    figure();% shows how much energy is contained in the null space of H'
    for j = 1:length(Yproj(:, 1))
        semilogx(f, 10*log10(P_nsm_stored(:, j)));
        hold on;
        semilogx(f, 10*log10(P_res_stored(:, j)));
    end
    xlabel('Frequency [Hz]');
    ylabel('[dB]');
    title('NSM Gravity Signal PSD');
    legend('projected signal', 'attitude error residual')
    grid on;
    
    figure()
    semilogx(f, 10*log10(E_ratio));
    hold on;
    xlabel('Frequency [Hz]');
    ylabel('NSM Transfer Function [dB]');
    title('Energy ratio NSM vs residual');
    grid on;

end

