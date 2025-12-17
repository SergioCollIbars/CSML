function [Q] = processNoise(QT, Qb, At, Ns)
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
   

    % process noise
    Q = zeros(Ns, Ns);

    % update Q
    I = eye(3, 3);
    Gamma = At * [At/2*I;I];
    q = Gamma * QT * Gamma';
    Nq = length(q);
    q = At.*[1/3*At^2*QT, 1/2*At*QT;...
        1/2*At*QT, QT];
    Q(1:Nq, 1:Nq) = q;

    Q(7:12, 1:12) = [zeros(6, 6), Qb];
end

