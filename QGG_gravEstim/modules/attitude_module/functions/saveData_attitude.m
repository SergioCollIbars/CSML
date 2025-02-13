function saveData_attitude(AttitudeObj, name)
    %%                       SAVE DATA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/11/2022
    %
    %   Description: This function save object data into txt file.
    %
    %   Input:
    %       AttitudeObj: attitude object
    %       name: file name
    %
    %   Output: null
    %
    % --------------------------------------------------------------------%

    % obtain complete file path
    path = "./data_files/" + name;

    % obtain variables to save
    t = AttitudeObj.t';
    
    theta1 = AttitudeObj.theta1';
    theta1Dot = AttitudeObj.theta1Dot';
    theta1Ddot = AttitudeObj.theta1Ddot';

    theta2 = AttitudeObj.theta2';
    theta2Dot = AttitudeObj.theta2Dot';
    theta2Ddot = AttitudeObj.theta2Ddot';

    theta3 = AttitudeObj.theta3';
    theta3Dot = AttitudeObj.theta3Dot';
    theta3Ddot = AttitudeObj.theta3Ddot';

    omega_x = AttitudeObj.omega(1, :)';
    omega_y = AttitudeObj.omega(2, :)';
    omega_z = AttitudeObj.omega(3, :)';

    Omega_x = AttitudeObj.Omega(1, :)';
    Omega_y = AttitudeObj.Omega(2, :)';
    Omega_z = AttitudeObj.Omega(3, :)';

    % create table
    T = table(t, theta1, theta1Dot, theta1Ddot, theta2, theta2Dot, ...
        theta2Ddot, theta3, theta3Dot, theta3Ddot, omega_x, omega_y, ...
        omega_z, Omega_x, Omega_y, Omega_z, 'VariableNames', ...
        {'TIME', 'theta1', 'theta1Dot', 'theta1Ddot', 'theta2', ...
        'theta2Dot', 'theta2Ddot', 'theta3', 'theta3Dot', 'theta3Ddot' ...
        'omega_x', 'omega_y', 'omega_z', 'Omega_x', 'Omega_y', 'Omega_z'});

    % save table
    writetable(T, path);
end