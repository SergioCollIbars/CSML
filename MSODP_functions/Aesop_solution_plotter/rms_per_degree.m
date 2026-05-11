function rms_deg = rms_per_degree(C_err, S_err)

    n_max = size(C_err,1) - 1;
    rms_deg = NaN(n_max+1,1);

    for n = 0:n_max
        cvals = C_err(n+1,1:n+1);
        svals = S_err(n+1,1:n+1);

        vals = [];

        % C coefficients: m = 0,...,n
        vals = [vals, cvals(~isnan(cvals))];

        % S coefficients: usually m = 1,...,n, since S_n0 = 0 / undefined
        if n >= 1
            s_use = svals(2:n+1);
            vals = [vals, s_use(~isnan(s_use))];
        end

        if ~isempty(vals)
            rms_deg(n+1) = sqrt(mean(vals.^2));
        end
    end
end