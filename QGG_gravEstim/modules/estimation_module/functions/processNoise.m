function [Q] = processNoise(QT, At, varObj)
%%                          PROCESS NOISE
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 03/05/2023
    %
    %   Description: Compute process noise time update based on the data
    %   gap. Assuming constant QTilde value.
    %
    %   Input: 
    %       QT: QT at the current time step
    %       DG: data gap time
    %       At: time span
    %
    %   Output:
    %       Q: process noise from ti 2 ti-1
    % --------------------------------------------------------------------%
    Q = zeros(varObj.Np, varObj.Np);

    % update Q
    I = eye(3, 3);
    Gamma = At * [At/2*I;I];
    q = Gamma * QT * Gamma';
    
    Nq = length(q);

    Q(1:Nq, 1:Nq) = q;
end

