function [] = plot_trajectory(time,state)
    % Plot the spacecraft trajectory.
    
    r = state(:, 1:3)'; v = state(:, 4:6)';
    
    % plot orbit radius & 3D trajectory
    figure()
    subplot(1, 2, 1)
    plot(time./3600, vecnorm(r)./1E3, 'LineWidth', 2)
    grid on;
    title('orbit radius norm')
    xlabel('Time [hr]')
    ylabel('[km]')

    subplot(1, 2, 2)
    plot3(r(1, :)./1E3, r(2, :)./1E3, r(3, :)./1E3, 'LineWidth', 2)
    grid on;
    title('3D trajectory ACI frame')
    xlabel('X [Km]')
    ylabel('Y [km]')
    zlabel('Z [Km]')
    hold all;

     % create planet surface
    [~,~,~] = sphere;
    scale =  450./1E3;  % Bennu object scale factor
    obj = readObj('Bennu-Radar.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 

    trisurf(f,p(:,1),p(:,2),p(:,3),'LineWidth', 0.1);
    colormap(gray);

    % plot axis
    mAxis = max(max(r./1E3));
    axis([0 mAxis 0 mAxis 0 mAxis])
    hold all;
    quiver3(0,0,-max(0),0,0,scale*(1.3),'r','LineWidth',1.5)
    quiver3(0,-max(0),0,0,scale*(1.3),0,'r','LineWidth',1.5)
    quiver3(-max(0),0,0,scale*(1.3),0,0,'r','LineWidth',1.5)
    text(0,0,scale*(1.3),'K','Color','r')
    text(0,scale*(1.3),0,'J','Color','r')
    text(scale*(1.3),0,0,'I','Color','r')

    axis equal
end

