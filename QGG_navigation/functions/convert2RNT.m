function [r_RTN,v_RTN, cov_RTN] = convert2RNT(r, v, t, P, Ns)
    r_RTN  = r .* 0;
    v_RTN  = v .* 0;
    cov_RTN = P .* 0;

    for j = 1:length(t)
        r_t = r(:, j); % inertial
        v_t = v(:, j); % inertial

        % compute rotation matrix
        [NB] = RTN2ECI(r_t, v_t);
        BN = NB';

        % convert states
        r_RTN(:, j) = BN * r_t;
        v_RTN(:, j) = BN * v_t;

        % convert uncertainty;
        if(Ns == 6)
            rotMat = [BN, zeros(3, 3);zeros(3, 3), BN];
        else
            rotMat = [BN, zeros(3, 4);zeros(3, 3), BN, zeros(3, 1); zeros(1, 7)];
        end
        p = reshape(P(j, 1:Ns*Ns), [Ns, Ns]);
        pp =  rotMat * p * rotMat';
        cov_RTN(j, 1:Ns*Ns) = reshape(pp, [1, Ns*Ns]);
    end
end

