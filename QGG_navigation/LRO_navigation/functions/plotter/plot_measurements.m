function plot_measurements(time, Y_N, bias, signal_error, ...
                           instrumentParams, Xf, orientation)
    % Plot gradiometer measurements and bias 

    % -------------------------
    % Compute measurements (body frame)
    % -------------------------
    BN_mat = compute_orientation_SC(time, Xf(1:6,:)', orientation);

    Y = compute_orientation_Meas(time, BN_mat, Y_N, signal_error);

    % Measurement mask
    mask = instrumentParams(:,1);
    gray = [.7 .7 .7];

    Nm = size(Y,1);
    tt = ["xx","xy","xz","yy","yz","zz"];

    utc  = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

    % -------------------------
    % UI figure + layout
    % -------------------------
    fig = uifigure('Name','Gradiometer Measurements','Position',[100 100 1200 800]);

    root = uigridlayout(fig,[1 1]);
    root.RowHeight    = {'1x'};
    root.ColumnWidth = {'1x'};
    root.Padding      = [0 0 0 0];

    tg = uitabgroup(root);

    tabMeas = uitab(tg,'Title','Measurements');
    tabBias = uitab(tg,'Title','Bias');

    % =========================================================
    % Tab 1: Measurements (2x3 grid)
    % =========================================================
    gl1 = uigridlayout(tabMeas,[2 3]);
    gl1.RowHeight    = {'1x','1x'};
    gl1.ColumnWidth = {'1x','1x','1x'};
    gl1.Padding      = [10 10 10 10];
    gl1.RowSpacing   = 10;
    gl1.ColumnSpacing = 10;

    for k = 1:Nm
        ax = uiaxes(gl1);

        if mask(k)==1
            color = 'b';
        else
            color = gray;
        end

        plot(ax, tUTC, Y(k,:)./1E3, 'LineWidth', 2, 'Color', color);
        grid(ax,'on');
        title(ax, tt(k));
        ylabel(ax,'[Eotvos]');
        xlabel(ax,'Epoch');
    end

    % =========================================================
    % Tab 2: Bias (2x3 grid)
    % =========================================================
    gl2 = uigridlayout(tabBias,[2 3]);
    gl2.RowHeight    = {'1x','1x'};
    gl2.ColumnWidth = {'1x','1x','1x'};
    gl2.Padding      = [10 10 10 10];
    gl2.RowSpacing   = 10;
    gl2.ColumnSpacing = 10;

    for k = 1:Nm
        ax = uiaxes(gl2);

        if mask(k)==1
            color = 'b';
        else
            color = gray;
        end

        plot(ax, tUTC, bias(k,:), 'LineWidth', 2, 'Color', color);
        grid(ax,'on');
        title(ax, tt(k));
        ylabel(ax,'[milli-Eotvos]');
        xlabel(ax,'Epoch');
    end

end

