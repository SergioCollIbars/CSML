%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class Definition: MultiBodySystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef MultiBodySystem
    properties
        % -----------------------------------------------------------------
        % System and normalization properties.
        % -----------------------------------------------------------------
        mu
        lStar
        tStar
        mStar
        
        % -----------------------------------------------------------------
        % Radius of each primary body and the location of each equilibrium
        % point.
        % -----------------------------------------------------------------
        p1R
        p2R
        posL1
        posL2
        posL3
        posL4
        posL5
    end
    
    methods
        % -----------------------------------------------------------------
        % Class Constructor
        % -----------------------------------------------------------------
        function obj = MultiBodySystem(systemName)
            params = systemParams(systemName);
            obj.mu = params(1);
            obj.lStar = params(2);
            obj.tStar = params(3);
            obj.mStar = params(4);
            obj.p1R = params(5);
            obj.p2R = params(6);
            
            obj.posL1 = computeLPos(1, obj.mu, 1e-15);
            obj.posL2 = computeLPos(2, obj.mu, 1e-15);
            obj.posL3 = computeLPos(3, obj.mu, 1e-15);
            obj.posL4 = computeLPos(4, obj.mu, 1e-15);
            obj.posL5 = computeLPos(5, obj.mu, 1e-15);
        end
    end
end