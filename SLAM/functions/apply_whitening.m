function [pref_W, Hc_W, H_HOT] = apply_whitening(R, pref_vec,...
    H_vec, Hn_vec, Hn2_vec)
    Nm = length(pref_vec(:, 1, 1));
    Nt = length(pref_vec(1, :));
    
    meas = pref_vec.*0;
    hc    = H_vec.*0; 
    hn    = Hn_vec.*0; 
    hn2   = Hn2_vec.*0;
    Ncs = length(H_vec(1,:, 1));
    for k = 1:Nm % weight measurements
        L = chol(R, 'lower');
        
        % residual withening
        meas(k, :) = (L \ pref_vec(k, :)')';
       
        % gravity coefficients withening
        H           = squeeze(H_vec(k, :, :));
        hc(k, :, :) =  (L \ H.').';

        % nuisance withening
        if(~isempty(Hn_vec))
            H           = squeeze(Hn_vec(k, :, :));
            hn(k, :, :) =  (L \ H.').';

            H            = squeeze(Hn2_vec(k, :, :));
            hn2(k, :, :) =  (L \ H.').';
        end 
    end

    % re-arrange gravity coefficient partials
    Hp = permute(hc, [1 3 2]);   % Now size = 6 x Nt x 45
    Hc_W = reshape(Hp, 6*Nt, Ncs);
    
    if(~isempty(Hn_vec))
        % re-arrange nuisance partials
        Hp = permute(hn, [1 3 2]);   % Now size = 6 x Nt x 3
        Hn_W = reshape(Hp, 6*Nt, 3);

        Hp = permute(hn2, [1 3 2]);   % Now size = 6 x Nt x 6
        Hn2_W = reshape(Hp, 6*Nt, 6);
    end

    % re-arrange measurement residuals
    pref_W         = reshape(meas, 6*Nt, 1);

    % eliminate nuisance parameters
    pref_W_new = nan(3*Nt, 1);
    Hc_W_new   = nan(3*Nt, Ncs); Hn_HOT_W_new = nan(3*Nt, 6);
    if(~isempty(Hn_vec))
        % Do block stacking
        Hn_block  = zeros(6*Nt, 3*Nt); 
        Hn2_block = zeros(6*Nt, 6*Nt); 
        for k = 1:Nt
% %             maxInd = 3 * k; minInd = maxInd - 2;
            maxInd2 = 6*k; minInd2 = maxInd2 - 5;
% %             Hn = Hn_W(minInd2:maxInd2, :);

            r = (k-1)*6 + (1:6);    % row indices for block k
            c = (k-1)*3 + (1:3);    % column indices for block k
            Hn_block(r, c)  = Hn_W(minInd2:maxInd2, :);

            r = (k-1)*6 + (1:6);    % row indices for block k
            c = (k-1)*6 + (1:6);    % column indices for block k
            Hn2_block(r,c) = Hn2_W(minInd2:maxInd2, :);
    
% %             C = null(Hn');
% %         
% %             pref_W_new(minInd:maxInd, 1)   = C' * pref_W(minInd2:maxInd2, 1);
% %             Hc_W_new(minInd:maxInd, :)     = C' * Hc_W(minInd2:maxInd2, :);
% % 
% %             Hn_HOT_W_new(minInd:maxInd, :) = C' * Hn2_W(minInd2:maxInd2, :);
        end
        
        % update residuals & partials without nuisance parameters effects
% %         pref_W = pref_W_new;
% %         Hc_W   = Hc_W_new;
% %         H_HOT  = Hn_HOT_W_new;
        
        C = null(Hn_block');

        pref_W = C' * pref_W;
        Hc_W   = C' * Hc_W;
        H_HOT  = C' * Hn2_block;
    end
end

