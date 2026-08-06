function plot_results(time, state_true, bias_true, SF_true, X_EKF, P_EKF, ...
                      posfit, mask, orientation, folder_name, NB_MOON)
    % Plot EKF results

    %% --- Time / sizes
    Nt = length(time);
    Nx = size(X_EKF,1);

    %% --- Orientation
    BN_mat = compute_orientation_SC(time, X_EKF(1:6,:)', orientation, ...
        NB_MOON);

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

    errB  = bias_true - X_EKF(7:12,:);
    errSF = SF_true.*0;

    %% --- UI figure + tabs
    fig_title = folder_name + '/Filter Results';
    fig = uifigure('Name',fig_title,'Position',[100 100 1200 800]);

    root = uigridlayout(fig,[1 1]);
    root.RowHeight = {'1x'};
    root.ColumnWidth = {'1x'};
    root.Padding = [0 0 0 0];

    tg = uitabgroup(root);

    tabP  = uitab(tg,'Title','Position Error');
    tabV  = uitab(tg,'Title','Velocity Error');
    tabB  = uitab(tg,'Title','Bias Error');
    tabS  = uitab(tg,'Title','SF Error');
    tabPo = uitab(tg,'Title','Posfit Measurement');

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
        
        rmsVal = rms(3*sigma(k,10:end), 'omitnan');

        plot(ax, tUTC, errorP(k,:), 'LineWidth', 2);
        plot(ax, tUTC, +3*sigma(k,:), 'k', 'LineWidth', 2);
        plot(ax, tUTC, -3*sigma(k,:), 'k', 'LineWidth', 2);
        ylim(ax, [-2*rmsVal 2*rmsVal]);

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
        
        rmsVal = rms(3*sigma(k+3,10:end).*1E3, 'omitnan');

        plot(ax, tUTC, errorV(k,:).*1E3, 'LineWidth', 2);
        plot(ax, tUTC, +3*sigma(k+3,:).*1E3, 'k', 'LineWidth', 2);
        plot(ax, tUTC, -3*sigma(k+3,:).*1E3, 'k', 'LineWidth', 2);
        ylim(ax, [-2*rmsVal 2*rmsVal]);

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

    %% =========================
    % Tab 4: SF error (2x3)
    %% =========================
    glS = uigridlayout(tabS,[2 3]);
    glS.RowHeight = {'1x','1x'};
    glS.ColumnWidth = {'1x','1x','1x'};
    glS.Padding = [10 10 10 10];
    glS.RowSpacing = 10;
    glS.ColumnSpacing = 10;

    ttB = ["xx","xy","xz","yy","yz","zz"];
    for k = 1:6
        ax = uiaxes(glS);
        hold(ax,'on');

        if mask(k)==1
            rmsVal = rms(3*sigma(12+k,10:end), 'omitnan');
            plot(ax, tUTC, errSF(k,:), 'LineWidth', 2);
            plot(ax, tUTC, +3*sigma(12+k,:), 'k', 'LineWidth', 2);
            plot(ax, tUTC, -3*sigma(12+k,:), 'k', 'LineWidth', 2);
            ylim(ax, [-2*rmsVal, 2*rmsVal]);
        else
            title(ax, ttB(k) + " (off)");
            grid(ax,'on');
            ax.XColor = [0.6 0.6 0.6];
            ax.YColor = [0.6 0.6 0.6];
            continue
        end

        grid(ax,'on');
        title(ax, ttB(k));
        ylabel(ax,'[-]');
        xlabel(ax,'Epoch');
    end

    %% =========================
    % Tab 5: Posit measurement (2x3) 
    %% =========================
    glPo = uigridlayout(tabPo,[2 3]);
    glPo.RowHeight = {'1x','1x'};
    glPo.ColumnWidth = {'1x','1x','1x'};
    glPo.Padding = [10 10 10 10];
    glPo.RowSpacing = 10;
    glPo.ColumnSpacing = 10;

    ttPo = ["xx","xy","xz","yy","yz","zz"];
    maskPositIdx = 1:6;  % <-- adjust if posit uses different mask indices

    for k = 1:6
        ax = uiaxes(glPo);
        hold(ax,'on');

        isActive = (k <= numel(maskPositIdx)) && (mask(k) == 1);

        if  isActive
            rmsVal = rms(posfit(k,:), 'omitnan');
            lgdStr = sprintf('RMS = %.3g',rmsVal);
            plot(ax, tUTC, posfit(k,:), 'LineWidth', 2, ...
                'DisplayName', lgdStr);
            legend(ax,'show','Location','best');
        else
            title(ax, ttPo(k) + " (off)");
            grid(ax,'on');
            ax.XColor = [0.6 0.6 0.6];
            ax.YColor = [0.6 0.6 0.6];
            continue
        end

        grid(ax,'on');
        title(ax, ttPo(k));
        ylabel(ax,'[mE]');
        xlabel(ax,'Epoch');
    end

end
