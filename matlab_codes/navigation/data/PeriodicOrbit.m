%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class Definition: PeriodicOrbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef PeriodicOrbit < handle
    properties
        % -----------------------------------------------------------------
        % System and orbit properties.
        % -----------------------------------------------------------------
        systemName
        mu
        centerIdx
        iState
        Xf
        period
        %% 
        CJ
        sIdx
        sVals
        sVecs
        monoMat
        timeC
        
        % -----------------------------------------------------------------
        % Discretized node properties.
        % -----------------------------------------------------------------
        sNodes
        tNodes
        eNodes
    end
    
    methods
        % -----------------------------------------------------------------
        % Class Constructor
        % -----------------------------------------------------------------
        function obj = PeriodicOrbit(systemName, centerIdx, orbit)
            obj.systemName = systemName;
            params = systemParams(obj.systemName);
            obj.mu = params(1);
            
            obj.centerIdx = centerIdx;
            obj.iState = orbit.iState;
            obj.Xf = orbit.Xf;
            obj.period = orbit.period;
            obj.CJ = orbit.CJ;
            obj.sIdx = orbit.sIdx;
            obj.sVals = orbit.sVals;
            obj.sVecs = orbit.sVecs;
            obj.monoMat = orbit.monoMat;
            obj.timeC = 1 / abs(real(log(max(abs(obj.sVals)))));
        end
        
        % -----------------------------------------------------------------
        % Output the mu value associated with the orbit's multi-body
        % system.
        % -----------------------------------------------------------------
        function mu = getMu(obj)
            mu = obj.mu;
        end
        
        % -----------------------------------------------------------------
        % Output the number of nodes.
        % -----------------------------------------------------------------
        function numNodes = getNumNodes(obj)
            numNodes = length(obj.sNodes(:,1));
        end
        
        % -----------------------------------------------------------------
        % Output the states, times, and event labels for all nodes.
        % -----------------------------------------------------------------
        function [sNodes, tNodes, eNodes] = getAllNodes(obj)
            sNodes = obj.sNodes;
            tNodes = obj.tNodes;
            eNodes = obj.eNodes;
        end
        
        % -----------------------------------------------------------------
        % Output the state, time, and event label for a desired node.
        % -----------------------------------------------------------------
        function [sNode, tNode, eNode] = getNode(obj, nodeIdx)
            if (nodeIdx > length(obj.sNodes(:,1)) || nodeIdx < 1)
                error("Invalid nodeIdx!");
            else
                sNode = obj.sNodes(nodeIdx,:);
                tNode = obj.tNodes(nodeIdx);
                eNode = obj.eNodes(nodeIdx);
            end
        end
        
        % -----------------------------------------------------------------
        % Output a desired state component for all nodes.
        % -----------------------------------------------------------------
        function sNodeComp = getNodeComp(obj, compIdx)
            if (compIdx > length(obj.sNodes(1,:)) || compIdx < 1)
                error("Invalid compIdx!");
            else
                sNodeComp = obj.sNodes(:,compIdx);
            end
        end
        
        % -----------------------------------------------------------------
        % Output the states of the nodes in a row vector.
        % -----------------------------------------------------------------
        function rowNodes = reshapeNodesRow(obj)
            rowNodes = reshape(obj.sNodes.', 1, length(obj.sNodes(:,1)) * length(obj.sNodes(1,:)));
        end
        
        % -----------------------------------------------------------------
        % Output the states of the nodes in a column vector.
        % -----------------------------------------------------------------
        function colNodes = reshapeNodesCol(obj)
            colNodes = reshape(obj.sNodes.', length(obj.sNodes(:,1)) * length(obj.sNodes(1,:)), 1);
        end
        
        % -----------------------------------------------------------------
        % Set the discretization of the periodic orbit.
        % -----------------------------------------------------------------
        function obj = setNodeSpace(obj, nodeSpaceIdx, numIntNodes)
            [obj.sNodes, obj.tNodes, obj.eNodes] = genNodesPO(obj.mu, obj.centerIdx, obj.Xf.', obj.period, nodeSpaceIdx, numIntNodes);
        end
        
        % -----------------------------------------------------------------
        % Generate an integrated trajectory from the current nodes.
        % -----------------------------------------------------------------
        function [t,q] = genTrajFromNodes(obj)
            t = [];
            q = [];
            numNodes = length(obj.tNodes);
            sc.mu = obj.mu;
            sc.stateI = obj.iState;
            for i = 1:(numNodes + 1)
                if (i == 1)
                    sc.timeSpan = [0; obj.tNodes(i)];
                elseif (i == (numNodes + 1))
                    sc.timeSpan = [obj.tNodes(i-1); obj.period];
                else
                    sc.timeSpan = [obj.tNodes(i-1); obj.tNodes(i)];
                end
                [tCurr, qCurr] = propCR3BP(sc, 0);
                sc.stateI = qCurr(end,:).';
                t = [t; tCurr];
                q = [q; qCurr];
            end
        end
    end
end