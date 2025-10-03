function saveData(t, state, X, meas)
    %%                    SAVE DATA FUNCTION
    % Description: save the computed estimation parameters in 
    % 'QGG_navigation.m' file.
    % Author: Sergio Coll Ibars
    % Date: 04/24/2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check vectors are in stack format [Nt, n];
    Nt = length(t);
    t = reshape(t, [Nt, 1]);
    
    posT = state(:, 1:3)';  % 'true' position state
    velT = state(:, 4:6)';  % 'true' velocity state

    posE = X(1:3, :);  % 'estimated' position state
    velE = X(4:6, :);  % 'estimated' velocity state

    % create and save data
    M = [t, posT', velT', posE', velE'];    % Dim: Nt x 24
    writematrix(M, "data/stateData.txt");

    M = [t, meas'];                         % Dim: Nt x 6
    writematrix(M, "data/measData.txt");
end

