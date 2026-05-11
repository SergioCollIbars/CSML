clc;
clear;
close all;
format long g;

%%                HARMONIC COEFICIENT PARSER FOR CLONES
% Description: Script to generate .txt files with the harmonic coefficents
% ordered in a single line using the SHA format reader provided by NASA.

% INPUT
addpath('/Users/sergiocollibars/Desktop/CSML/matlab_codes/QGG_equations/functions');
out_path   = 'data_out/';
input_path = 'data_in/';

files  = dir(input_path);
files = files(~[files.isdir]);    % keep only files (remove folders)

R = 1738E3;                        % Reference Radius [m]
for f = 1:length(files)
    disp('Processing: '); disp( string(files(f).name));
    [~,~,ext] = fileparts(files(f).name);

    isTxt = strcmpi(ext, '.txt');
    if(~isTxt)
        continue
    end

    file_1200  = readmatrix(input_path + string(files(f).name));
    
    % count coeff
    [Nc_1200, Ns_1200]       = countCoeff(1200);
    
    Nc_col    = 3;     % line to extract C coeff
    Ns_col    = 4;     % line to extract S coeff
    init_line = 1;     % line to start file reading
    [M] = orderMatrix(file_1200, Nc_1200, Ns_1200, Nc_col, Ns_col,...
        init_line, 0);
    
    % save matrix into txt file
    T = table([1200; R; 1; M],...
        'VariableNames', {'VALUES'});
    
    % save table
    writetable(T, out_path +  string(files(f).name));
end
disp('FINISHED!')

%% FUNCTIONS

function [Nc, Ns] = countCoeff(n)
    % count number of C and S coeff
    Nc = 1;
    for k = 2:n
        Nc = Nc + k + 1;
    end
    Ns = 0;
    for k = 2:n
        Ns = Ns + k;
    end
end

function [N] = NormFactor(n, m)
    % Description: given degree, n and order, m compute the normalice
    % factor
    if(m == 0)
        delta = 1;
    else
        delta = 0;
    end
    fac1 = factorial(n - m);
    fac2 = factorial(n + m);
    N = ((2 - delta)*(2*n + 1) * fac1 /fac2)^(0.5);
end

function [Matrix] = orderMatrix(file, Nc, Ns, Nc_col, Ns_col,...
    init_line, normalized)
    % Order C and S matrices in a 1x1 matrix
    Matrix = ones(Nc+Ns, 1) * NaN;

    % Point mass solution
    Matrix(1, 1) = 0;

    % extract values
    Cval = file(init_line:end, Nc_col);
    Sval = file(init_line:end, Ns_col);
    
    n = 2;
    m = 0;
    countFile = 1;
    for j = 2:Nc
        if(normalized == 1)
            Norm = NormFactor(n, m);
        else
            Norm = 1;           % fully normalized
            %Norm = sqrt(4*pi);  % fully 4pi normalized
        end

        Matrix(j) = Cval(countFile) * Norm;
        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n +1;
        end
        countFile = countFile + 1;
    end

    n = 2;
    m = 0;
    count = Nc + 1;
    for j = 1:Nc
        if(normalized == 1)
            Norm = NormFactor(n, m);
        else
            Norm = 1;
        end
        
        if(m ~= 0)
            Matrix(count) = Sval(j) * Norm;
            count = count + 1;
        end

        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
    end
end