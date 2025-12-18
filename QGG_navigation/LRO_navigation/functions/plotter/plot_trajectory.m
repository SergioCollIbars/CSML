function plot_trajectory(time, state)

    % Data (unchanged)
    r = state(:, 1:3)'; 
    v = state(:, 4:6)';

    utc  = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

    R   = cspice_bodvrd('MOON', 'RADII', 3) * 1E3;
    R_M = R(1);

    H = vecnorm(r) - R_M;

    % ---- UI figure ----
    fig = uifigure('Name','Trajectory Plots','Position',[100 100 1200 800]);

    % Make a layout that fills the whole window
    root = uigridlayout(fig,[1 1]);
    root.RowHeight = {'1x'};
    root.ColumnWidth = {'1x'};
    root.Padding = [0 0 0 0];

    % Put the tab group inside the root layout (so it fills)
    tg = uitabgroup(root);

    tab1 = uitab(tg,'Title','Time Series');
    tab2 = uitab(tg,'Title','3D Trajectory');
    tab3 = uitab(tg,'Title','Alt');

    %% --- Tab 1
    gl1 = uigridlayout(tab1,[1 2]);
    gl1.RowHeight = {'1x'};
    gl1.ColumnWidth = {'1x','1x'};
    gl1.Padding = [10 10 10 10];
    gl1.ColumnSpacing = 10;

    ax1 = uiaxes(gl1);
    ax2 = uiaxes(gl1);

    plot(ax1, tUTC, vecnorm(r)./1E3, 'LineWidth', 2);
    grid(ax1,'on'); title(ax1,'orbit radius norm');
    xlabel(ax1,'Epoch'); ylabel(ax1,'[km]');

    plot(ax2, tUTC, vecnorm(v), 'LineWidth', 2);
    grid(ax2,'on'); title(ax2,'orbit velocity norm');
    xlabel(ax2,'Epoch'); ylabel(ax2,'[m/s]');

    %% --- Tab 2
    gl2 = uigridlayout(tab2,[1 1]);
    gl2.RowHeight = {'1x'};
    gl2.ColumnWidth = {'1x'};
    gl2.Padding = [10 10 10 10];

    ax3 = uiaxes(gl2);
    hold(ax3,'on');

    plot3(ax3, r(1,:), r(2,:), r(3,:), 'LineWidth', 2);

    [Xs,Ys,Zs] = sphere(100);
    surf(ax3, R_M*Xs, R_M*Ys, R_M*Zs, ...
        'EdgeColor','none', 'FaceColor',[0.8 0.8 0.8]);

    axis(ax3,'equal'); grid(ax3,'on');
    xlabel(ax3,'X [m]'); ylabel(ax3,'Y [m]'); zlabel(ax3,'Z [m]');
    title(ax3,'LRO Trajectory around the Moon');

    %% --- Tab 3
    gl3 = uigridlayout(tab3,[1 1]);
    gl3.RowHeight = {'1x'};
    gl3.ColumnWidth = {'1x'};
    gl3.Padding = [10 10 10 10];

    ax4 = uiaxes(gl3);

    plot(ax4, tUTC, H./1E3, 'LineWidth', 2);
    grid(ax4,'on');
    xlabel(ax4,'Epoch'); ylabel(ax4,'[km]');
    title(ax4,'orbit altitude');

end



