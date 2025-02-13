function [] = plot_trajectory(state_t, body)
    figure()
    plot3(state_t(:,1), state_t(:,2), state_t(:,3), 'LineWidth', 2, ...
        'Color','r')
    axis equal;
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title('S/C trajectory')
    grid on;
    
    % plot axis
    mAxis = max(max(state_t(:, 1:3)));
    axis([0 mAxis 0 mAxis 0 mAxis])
    hold all;
    quiver3(0,0,-max(0),0,0,max(zlim),'b','LineWidth',1)
    quiver3(0,-max(0),0,0,max(ylim),0,'b','LineWidth',1)
    quiver3(-max(0),0,0,max(xlim),0,0,'b','LineWidth',1)
    text(0,0,max(zlim),'K','Color','b')
    text(0,max(ylim),0,'J','Color','b')
    text(max(xlim),0,0,'I','Color','b')
    
    scale =  450;  % Bennu object scale factor
    if(body == "BENNU")
        obj = readObj('Bennu-Radar.obj');
    elseif( body == "EROS")
        obj = readObj('Eros-Poly.obj');
    end
    p = obj.v * 2 * scale;
    f = obj.f.v ; 
    
    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
end
