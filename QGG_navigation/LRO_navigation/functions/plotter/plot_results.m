function plot_results(time, state_true, bias_true, X_EKF, P_EKF, ...
                      mask, orientation)
    % Plot EKF results 

    %% --- Time / sizes
    Nt = length(time);
    Nx = size(X_EKF,1);

    %% --- Orientation
    BN_mat = compute_orientation_SC(time, X_EKF(1:6,:)', orientation);

    %% --- Convert time to UTC
    utc  = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

    %% --- Errors and sigmas
    sigma  = nan(Nx, Nt);

    errorP = state_true(:,1:3)' - X_EKF(1:3,:);
    errorV = state_true(:,4:6)' - X_EKF(4:6,:);

    for k = 1:Nt
        p = reshape(P_EKF(k,:), [Nx, Nx]);

        maxInd = 3*k; 
        minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);

        A   = blkdiag(BN, BN, eye(6));
        P_B = A * p * A';
        sigma(:,k) = sqrt(diag(P_B));

        errorP(:,k) = BN * errorP(:,k);
        errorV(:,k) = BN * errorV(:,k);
    end

    errB = bias_true - X_EKF(7:end,:);

    %% --- UI figure + tabs
    fig = uifigure('Name','Filter Results','Position',[100 100 1200 800]);

    root = uigridlayout(fig,[1 1]);
    root.RowHeight = {'1x'};
    root.ColumnWidth = {'1x'};
    root.Padding = [0 0 0 0];

    tg = uitabgroup(root);

    tabP = uitab(tg,'Title','Position Error');
    tabV = uitab(tg,'Title','Velocity Error');
    tabB = uitab(tg,'Title','Bias Error');

    %% =========================
    % Tab 1: Position error
    %% =========================
    glP = uigridlayout(tabP,[1 3]);
    glP.RowHeight = {'1x'};
    glP.ColumnWidth = {'1x','1x','1x'};
    glP.Padding = [10 10 10 10];
    glP.ColumnSpacing = 10;

    ttP = ["X","Y","Z"];
    for k = 1:3
        ax = uiaxes(glP);
        hold(ax,'on');

        plot(ax, tUTC, errorP(k,:), 'LineWidth', 2);
        plot(ax, tUTC, +3*sigma(k,:), 'k', 'LineWidth', 2);
        plot(ax, tUTC, -3*sigma(k,:), 'k', 'LineWidth', 2);

        grid(ax,'on');
        title(ax, ttP(k));
        ylabel(ax,'[m]');
        xlabel(ax,'Epoch');
    end

    %% =========================
    % Tab 2: Velocity error
    %% =========================
    glV = uigridlayout(tabV,[1 3]);
    glV.RowHeight = {'1x'};
    glV.ColumnWidth = {'1x','1x','1x'};
    glV.Padding = [10 10 10 10];
    glV.ColumnSpacing = 10;

    ttV = ["V_x","V_y","V_z"];
    for k = 1:3
        ax = uiaxes(glV);
        hold(ax,'on');

        plot(ax, tUTC, errorV(k,:).*1E3, 'LineWidth', 2);
        plot(ax, tUTC, +3*sigma(k+3,:).*1E3, 'k', 'LineWidth', 2);
        plot(ax, tUTC, -3*sigma(k+3,:).*1E3, 'k', 'LineWidth', 2);
        ylim(ax, [-1 1]);

        grid(ax,'on');
        title(ax, ttV(k));
        ylabel(ax,'[mm/s]');
        xlabel(ax,'Epoch');
    end

    %% =========================
    % Tab 3: Bias error (2x3)
    %% =========================
    glB = uigridlayout(tabB,[2 3]);
    glB.RowHeight = {'1x','1x'};
    glB.ColumnWidth = {'1x','1x','1x'};
    glB.Padding = [10 10 10 10];
    glB.RowSpacing = 10;
    glB.ColumnSpacing = 10;

    ttB = ["xx","xy","xz","yy","yz","zz"];
    for k = 1:6
        ax = uiaxes(glB);
        hold(ax,'on');

        if mask(k)==1
            plot(ax, tUTC, errB(k,:), 'LineWidth', 2);
            plot(ax, tUTC, +3*sigma(6+k,:), 'k', 'LineWidth', 2);
            plot(ax, tUTC, -3*sigma(6+k,:), 'k', 'LineWidth', 2);
        else
            % If masked off, keep subplot but make it visually "inactive"
            title(ax, ttB(k) + " (off)");
            grid(ax,'on');
            ax.XColor = [0.6 0.6 0.6];
            ax.YColor = [0.6 0.6 0.6];
            continue
        end

        grid(ax,'on');
        title(ax, ttB(k));
        ylabel(ax,'[mE]');
        xlabel(ax,'Epoch');
    end

end

