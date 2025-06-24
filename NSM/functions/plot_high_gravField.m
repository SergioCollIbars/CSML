function [] = plot_high_gravField(n_max, Nc, Ns, value, tt)
    % Description: Compute the value vector (Nc+Ns x 1) in a pyramide plot.
    % Author: Sergio Coll Ibars
    % Date: 06/23/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [A, B] = list2mat(n_max, Nc, Ns, value);
    A(1,1) = value(1);

    J = ones(n_max+1, (n_max+1)*2) * NaN; S = fliplr(B);
    J(:, 1:n_max+1) = S;
    J(:, n_max+2:end) = A;

    figure;
    cmap = turbo(256);  % Original colormap with 256 colors
    im = imagesc(J);
    clim([0 100]);  % Force color scaling between 1% and 100%
    colormap(cmap);
    set(im, 'AlphaData', ~isnan(J)); % Make NaNs transparent
    hold on;

    c = colorbar;
    c.Ticks = [0, 20, 40, 60, 80, 100]; 
    c.TickLabels = {'1', '20', '40', '60', '80', '100 %'};
   
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

