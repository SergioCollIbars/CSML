function [estimObj, maxIter] = getParams(estimObj, OrbitObj, normalized,...
    mode)
%%                   GET PARAMETERS FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 11/01/2023
    %
    %   Description: This function returns the main parameters for the
    %   specified problem: Earth or Juno
    %
    %   Input:
    %       estimObj: estimation class object
    %       normalized: 1 for normalized coefficients, 0 otherwise
    %       orbitObj: orbital class object
    %       mode: extra parameter mode
    %  
    %   Output: 
    %        estimObj: estimation class object
    % --------------------------------------------------------------------%
    
    % read max iterations
    maxIter = estimObj.EstimData(6);
    
    % load extra parameter mode
    estimObj.extraParam = mode;

    % GM parameter
    estimObj.GM = estimObj.EstimData(8);

    % planet radius
    estimObj.Re = OrbitObj.OrbitData(10);

    % define SH degree
    n = estimObj.EstimData(2);
    Nc  = -2;       % # number of Cnm parameters
    for j = 1:n+1
        Nc = Nc + j;
    end
    Ns = Nc - n;    % # number of Snm parameters
    estimObj.Nc = Nc;
    estimObj.Ns = Ns;
    estimObj.n_max = n;

    % define state vector
    X = zeros(Nc+Ns, 1);

    % define extra parameters vector
    [E, E0, ~, ~, ~] = load_extraParams(mode, Nc, Ns);

    % initialize states
    [Cnm0, Snm0] = getkaula(n, Nc, Ns, normalized);
% %     Cnm0 = load("C.mat").C_mat';
% %     Snm0 = load("S.mat").S_mat';
% %     Cnm0(1) = 1;
% %     
    X0 = [Cnm0; Snm0];
    if(any(E0))
        X = [X; E];
        X0 = [X0; E0];
    end

    % save variables in the estimation object
    estimObj.X0 = X0;
    estimObj.delta_X0 = zeros(length(X0), 1);

    % pole parameters vector
    W = OrbitObj.W;
    W0 = OrbitObj.W0;
    RA = OrbitObj.RA;
    DEC = OrbitObj.DEC;
    estimObj.poleParams = [W, W0, RA, DEC];
end

