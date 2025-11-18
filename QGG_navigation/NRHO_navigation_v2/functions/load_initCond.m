function [X0] = load_initCond(planetParams, TIME)
    %%                    LOAD INITIAL CONDITIONS FUNCTION
    % Description: Based on the especified system, load the initial
    % conditions of the orbit.
    % Author: Sergio Coll Ibars
    % Date: 03/27/2024

    et = TIME(1)/planetParams(3);   % [-]
    [state, ~] = cspice_spkezr('-60000', et, 'J2000', 'NONE', '3');

    r0 = state(1:3);
    v0 = state(4:6);

    r0 = r0*1E3/planetParams(2);
    v0 = v0*1E3/(planetParams(3)*planetParams(2));

    % define initial conditions
    X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
end

