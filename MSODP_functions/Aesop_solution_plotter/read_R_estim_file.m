function [C_true, C_err, S_true, S_err] = read_R_estim_file(filename, n_max, SCALE)

    fid = fopen(filename, 'r');
    txt = fread(fid, '*char')';
    fclose(fid);

    % Extract only floating-point numbers
    float_tokens = regexp(txt, '[-+]?(?:\d+\.\d*|\.\d+)(?:[eEdD][-+]?\d+)?', 'match');
    vals = str2double(float_tokens);
    vals = vals(~isnan(vals));

    % Each coefficient: true, err, est
    if mod(numel(vals),3) ~= 0
        error('Float count not multiple of 3');
    end

    A = reshape(vals,3,[])';
    true_vals = A(:,1);
    err_vals  = A(:,2);

    % Expected number of coefficients
    N_expected = (n_max - 1); % zonals
    for n = 2:n_max
        N_expected = N_expected + 2*n; % m=1:n → C and S
    end

    % Trim if needed
    true_vals = true_vals(1:N_expected);
    err_vals  = err_vals(1:N_expected);

    % Allocate
    C_true = NaN(n_max+1, n_max+1);
    C_err  = NaN(n_max+1, n_max+1);
    S_true = NaN(n_max+1, n_max+1);
    S_err  = NaN(n_max+1, n_max+1);

    k = 1;

    % -----------------------
    % 1) Zonal: C_n0
    % -----------------------
    for n = 2:n_max
        C_true(n+1,1) = true_vals(k);
        C_err(n+1,1)  = err_vals(k);
        k = k + 1;
    end

    % -----------------------
    % 2) Degree-wise loop
    % -----------------------
    for n = 2:n_max
        for m = 1:n
            % C_nm
            C_true(n+1,m+1) = true_vals(k);
            C_err(n+1,m+1)  = err_vals(k);
            k = k + 1;

            % S_nm
            S_true(n+1,m+1) = true_vals(k);
            S_err(n+1,m+1)  = err_vals(k);
            k = k + 1;
        end
    end

    % Apply SCALE factor 
    C_true = C_true.*SCALE;
    S_true = S_true.*SCALE;

    C_err = C_err.*SCALE;
    S_err = S_err.*SCALE;
end