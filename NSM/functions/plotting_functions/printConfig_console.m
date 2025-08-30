function [] = printConfig_console(Planet, Solver, Errors, measMode, ...
    saveData, type_Attitude, type_posErr, type_attErr)
    % Print input configuration
    fprintf('\nSimulation Input Parameters:\n');
    fprintf('-----------------------------\n');
    fprintf('%-10s : %s\n', 'Planet',   Planet);
    fprintf('%-10s : %s\n', 'Solver',   Solver);
    fprintf('%-10s : %s\n', 'Errors',   Errors);
    if(measMode == "1")
        measMode = "GG + GYRO + ST";
        fprintf('%-10s : %s\n', 'measMode', measMode);
    else
        measMode = "GG + ST";
        fprintf('%-10s : %s\n', 'measMode', measMode);
    end
    if(saveData)
        name = "ON";
        fprintf('%-10s : %s\n', 'saveData', name);
    else
        name = "OFF";
        fprintf('%-10s : %s\n', 'saveData', name);
    end
    fprintf('%-10s : %s\n', 'Attitude orientation',   type_Attitude);
    fprintf('%-10s : %s\n', 'Attitude errors',   type_attErr);
    fprintf('%-10s : %s\n', 'Position errors',   type_posErr);
    fprintf('\n');
end