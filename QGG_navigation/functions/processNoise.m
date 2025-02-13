function [Q] = processNoise(QT, DG, At, Bw, type, Ns)
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
    
    % max DG allowed
    DG_max = 20; % [sec]

    % process noise
    Q = zeros(Ns, Ns);

    % update Q
    if(type == "SNC")
        if(DG < DG_max)
            I = eye(3, 3);
            Gamma = At * [At/2*I;I];
            q = Gamma * QT * Gamma';
            Nq = length(q);

            Q(1:Nq, 1:Nq) = q;
        end
    end
    
    if(type == "DMC")
        if(DG < DG_max)
            beta = Bw(1, 1);
            variance = QT;
            Q = DMCCovariance(At, variance, beta);
        end
    end
end

