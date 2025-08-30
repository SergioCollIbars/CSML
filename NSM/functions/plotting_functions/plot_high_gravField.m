function [] = plot_high_gravField(n_max, Nc, Ns, value, tt, limit, limitTicks, scale, digits)
    % Description: Compute the value vector (Nc+Ns x 1) in a pyramide plot.
    % Author: Sergio Coll Ibars
    % Date: 06/23/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [A, B] = list2mat(n_max, Nc, Ns, value);
    A(1,1) = value(1);

    J = ones(n_max+1, (n_max+1)*2) * NaN; S = fliplr(B);
    J(:, 1:n_max+1) = S;
    J(:, n_max+2:end) = A;
    J(J == 0) = NaN;    % set 0's to NaNs

%     figure;
    cmap = turbo(256);  % Original colormap with 256 colors
    im = imagesc(J);
    % % clim([1E-2 50]);      % Force color scaling between 1% and 100%
    clim([limit(1) limit(2)]);      % Force color scaling between 1% and 100%
    colormap(cmap);
    set(im, 'AlphaData', ~isnan(J)); % Make NaNs transparent
    hold on;
    set(gca, 'ColorScale', scale)

    c = colorbar;
    % % c.Ticks = [1E-2, 1E-1, 1, 5, 20, 50]; 
    c.Ticks = limitTicks; 
    c.TickLabels = ...
        arrayfun(@(x) num2str(x, digits), c.Ticks, 'UniformOutput', false);
   
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

