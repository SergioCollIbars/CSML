clc;
clear;
close all;
format long g;

%%                  HARMONIC COEFICIENT PARSER
% Description: Script to generate .txt files with the harmonic coefficents
% ordered in a single line using the SHA format reader provided by NASA.

% INPUT
addpath('./functions/');
addpath('//Users/sergiocollibars/Desktop/');
out_path   = './HARMCOEFS_MOON_1.txt';
input_path = 'GRGM660PRIM.txt';

normalized_ouput = 0;          % do you want normalized output coeff?
normalized       = 1;          % output file normalized?
n = 10;                        % max zonal harmonic
R = 1738E3;                    % Reference Radius [m]


% count coeff
[Nc, Ns] = countCoeff(n);

%%file = table2array(readtable(input_path));
file = readmatrix(input_path);

Nc_col    = 3;     % line to extract C coeff
Ns_col    = 4;     % line to extract S coeff
init_line = 3;     % line to start file reading
[M] = orderMatrix(file, Nc, Ns, Nc_col, Ns_col, init_line, ...
    normalized_ouput);


% save matrix into txt file
T = table([n; R; normalized; M],...
    'VariableNames', {'VALUES'});

% save table
writetable(T, out_path);

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
    Matrix(1, 1) = 1;

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