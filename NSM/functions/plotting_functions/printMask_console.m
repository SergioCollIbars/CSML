function [] = printMask_console(mask)
    % Measurement component labels
    labels = {'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'};
    
    % Print header
    fprintf('\nMeasurement Mask:\n');
    fprintf('------------------\n');
    fprintf('%4s', labels{:});
    fprintf('\n');
    
    % Print mask values
    fprintf('%4d', mask);
    fprintf('\n\n');
end
