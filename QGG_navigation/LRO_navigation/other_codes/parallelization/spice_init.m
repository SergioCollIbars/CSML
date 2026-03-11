function spice_init(metaKernelPath)
    % Avoid re-loading every time:
    persistent isInit
    if isempty(isInit) || ~isInit
        cspice_kclear;
        cspice_furnsh(metaKernelPath);
        isInit = true;
    end
end
