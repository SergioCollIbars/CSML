function [] = plot_signal(dataVec, t)
    dt = t(2) - t(1);
    fs = 1/dt;         % Hz

    % gradient color 
    cl = 'b';
    % plot gradiometer signal time series
    tt = ["\Gamma_{xx}", "\Gamma_{xy}", "\Gamma_{xz}", ...
        "\Gamma_{yy}", "\Gamma_{yz}", "\Gamma_{zz}"];

    % convert GPS time to UTC
    gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
    t_UTC = gps_epoch + seconds(t);        % date time 
    
    data = [dataVec(1, :);dataVec(2, :);dataVec(3,:);dataVec(5,:);...
        dataVec(6, :);dataVec(9, :)];
    figure()
    for j = 1:6
        subplot(2,3,j)

        [data(j, :)] = checkOutliers(data(j, :));

        plot(t_UTC, data(j, :)./1E-9, 'LineWidth', 2, ...
            'Color', cl);
        % % maxVal = max(data(j,:)./1E-9); minVal = min(data(j,:)./1E-9);
        % % ylim([2*minVal, 2*maxVal])
        grid on;
        ylabel('[Eotvos]')
        title(tt(j))
        %xlim([datetime(2012,11,02,0,0,0), datetime(2012,11,03,0,0,0)])
    end
    
% %     % Compute PSD using Welch's method
% %     for j = 1:6
% %         figure()
% %         % remove NaNs
% %         x = data(j, :)./1E-9;
% %         val = x(~isnan(x));
% %         [psd_vals, f] = pwelch(val);
% % 
% %         loglog(f, sqrt(psd_vals), 'LineWidth', 2); % sqrt to show in [E/√Hz]
% %         title(tt(j))
% %     end
end


% fucntion to check for outliers
function [dataOut] = checkOutliers(x)  
    Q1 = quantile(x, 0.25);
    Q3 = quantile(x, 0.75);
    IQR = Q3 - Q1;
    
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    
    outliers = (x < lower_bound) | (x > upper_bound);
    x(outliers) = NaN;
    dataOut = x;
end
