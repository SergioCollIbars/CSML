function [yvals, xvals] = orderValues(vals, n_max)
    [Nc, ~, Ncs] = count_num_coeff(n_max); 
    C_vals = vals(1:Nc-1);
    S_vals = vals(Nc:Ncs-1);

    % get values starting at n = 2
    yvals = []; xvals = [];
    for n = 2:n_max
        [Nc, Ns, ~] = count_num_coeff(n); 

        C_maxIdx = Nc-1;
        C_minIdx = (Nc-1) - n;

        S_maxIdx = Ns;
        S_minIdx = Ns - (n-1); 

        vec = [C_vals(C_minIdx:C_maxIdx); S_vals(S_minIdx:S_maxIdx)];
        yvals = [yvals; vec];

        Nvals = length(vec);
        Nmin = n;
        Nmax = n + 1;
        xvals    =  [xvals, linspace(Nmin, Nmax-0.1, Nvals)];
    end
end