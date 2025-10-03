function [T] = compute_meas_bias(meas, q, tau, TIME)
    %%                    COMPUTE MEASUREMENTS BIAS
    % Description: compute measurement bias and add them to the
    % measurements.
    % Author: Sergio Coll Ibars
    % Date: 04/09/2025
    % Data: meas: measurement vector
    %       q: noise PSD intensity
    %       tau: correlation parameter

    Nt = length(meas(1,:));
    Nm = length(meas(:, 1));
    dt = TIME(2) - TIME(1); % [-]

    % bias variance
    s = q * tau / 2 * (1 - exp(-2/tau * dt));

    % bias vector
    b = zeros(Nm, Nt);

    for j = 1:Nt-1
         w = normrnd(0, sqrt(s), [Nm, 1]);
         b(:, j + 1)    = exp(-dt/tau) * b(:, j) + w;      % [1 / s^2]
    end

    T = meas + b;
end

