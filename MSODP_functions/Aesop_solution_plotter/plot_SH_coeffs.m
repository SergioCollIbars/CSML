function [] = plot_SH_coeffs(C_nm,S_nm, n_max, tt, limits, scaleType)
    figure();
    A = abs(C_nm); B = abs(S_nm);
    J = ones(n_max+1, (n_max+1)*2) * NaN; S = fliplr(B);
    J(:, 1:n_max+1) = S;
    J(:, n_max+2:end) = A;
    J(J == 0) = NaN;    % set 0's to NaNs
    
    cmap = jet(256);               % Original colormap with 256 colors
    im = imagesc(J);
    if(~isempty(limits))
        clim(limits);
    end
    colormap(cmap);
    set(im, 'AlphaData', ~isnan(J)); % Make NaNs transparent
    hold on;
    set(gca, 'ColorScale', scaleType);
    
    c = colorbar;
    % % c.Ticks = [1E-2, 1E-1, 1, 5, 20, 50]; 
    % % c.Ticks = limitTicks; 
    % % c.TickLabels = ...
    % %     arrayfun(@(x) num2str(x, digits), c.Ticks, 'UniformOutput', false);
    
    % Set axis properties
    yticks(linspace(1, n_max + 1, 7));
    yticklabels(compose('%i', linspace(0, n_max, 7)));
    xticks(linspace(1, 2*n_max + 2, 7));
    xticklabels(compose('%i', linspace(-n_max, n_max, 7)));
    
    xlabel('Order');
    ylabel('Degree');
    hold off;
    title(tt);
end

