clear;
clc;
close all;
%%                  PLOT CORRELATIONS OVER TIME
% Date: 11/11/2025
% Author: Sergio Coll-Ibars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
data = load('cov_Mc_1.mat').data;
F_mat = data{1};
time  = data{2}; 

Nt_array = sum(~isnan(F_mat));     Nt = Nt_array(1);
Ns_array = sum(~isnan(F_mat), 2);  Ns = sqrt(Ns_array(1));

P = F_mat(1:Nt, 1:Ns*Ns);
fps = 120;
outFile = "covariance_time.mp4";

% create animation
makeCorrelationVideo(P, time, outFile, fps, Ns, Nt)


%% FUNCTIONS
function makeCorrelationVideo(P, t, outFile, fps, Ns, Nt)
% makeCorrelationVideo  Create a video of correlation matrices over time
%
%   makeCorrelationVideo(P, t, outFile, fps)
%
%   Inputs
%     P        : N x n t covariance matrices over time
%     t        : Nt-length time vector (datetime or numeric)
%     outFile  : output video filename (e.g., 'corr_video.mp4')
%     fps      : (optional) frames per second (default: 10)
%
%   Example
%     % P is n x n x Nt, t is 1xNt datetime vector
%     makeCorrelationVideo(P, t, 'corr.mp4', 12);

    if nargin < 4 || isempty(fps), fps = 10; end

    % Prepare video writer
    vw = VideoWriter(outFile, 'MPEG-4');
    vw.FrameRate = fps;
    vw.Quality   = 95;
    open(vw);

    % Figure setup (consistent size for all frames)
    fig = figure('Color','w','Units','pixels','Position',[100 100 800 700]);
    ax  = axes('Parent',fig);
    
    % Pre-create image for speed; we will update CData each frame
    hImg = imagesc(zeros(size(P,1), size(P,2)), 'Parent', ax);
    axis(ax,'image');
    set(ax,'YDir','normal');
    grid(ax,'on');
    ax.GridAlpha = 0.15;
    clim(ax,[0 1]);     % correlation range
    cb = colorbar(ax);
    cb.Label.String = 'Correlation';
    xlabel(ax,'State index');
    ylabel(ax,'State index');
    colormap(ax, turbo);

    % Title helper
    function txt = makeTitle(tt, k, Nt)
        if isdatetime(tt)
            txt = sprintf('Correlation  |  %s  |  frame %d/%d',...
                datestr(tt,'yyyy-mm-dd HH:MM:SS'), k, Nt);
        else
            txt = sprintf('Correlation  |  t = %.6g  |  frame %d/%d', tt, k, Nt);
        end
    end

    % Render frames
    for k = 1:Nt
        Ck = reshape(P(k, :), [Ns, Ns]);

        % Ensure symmetry (in case of numeric drift)
        Ck = 0.5*(Ck + Ck.');

        % Convert covariance to correlation robustly:
        d  = sqrt(max(diag(Ck), 0));         % avoid negative due to roundoff
        d(d==0) = eps;                       % avoid divide-by-zero
        denom = d * d.';                     % outer product
        Rk = abs(Ck ./ denom);                    % correlation matrix

        % Clip any tiny numeric spillover beyond [-1,1]
        Rk(Rk >  1) =  1;
        Rk(Rk < -1) = -1;

        % Update image
        set(hImg, 'CData', Rk);
        title(ax, makeTitle(t(k), k, Nt), 'FontWeight','bold');

        drawnow;
        frame = getframe(fig);
        writeVideo(vw, frame);
    end

    close(vw);
    close(fig);

    fprintf('Saved video: %s\n', outFile);
end
