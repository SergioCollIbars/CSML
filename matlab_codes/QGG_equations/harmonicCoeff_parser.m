clc;
clear;
close all;
format long g;

%%                  HARMONIC COEFICIENT PARSER
% Description: Script to generate .txt files with the harmonic coefficents
% ordered in a single line.

% INPUT
addpath('//Users/sergiocollibars/Desktop/CSML/codes/QGG/data_files/');
addpath('//Users/sergiocollibars/Desktop/Harmonics_Bennu_cd/');
addpath('//Users/sergiocollibars/Desktop/');
out_path = './HARMCOEFS_EROS_1.txt';

truth_values = "file";          % file / matrix
normalized = 1;                 % input file normalized?

% READ TRUTH VALUES
if(truth_values == "matrix")
    n = 6;
    R = 246;

    % count coeff
    [Nc, Ns] = countCoeff(n);

    C_Bennu = [1 0 0 0 0 0 0;...
                0 0 0 0 0 0 0;...
                -0.017511 -0.0000023 0.0058194 0 0 0 0;...
                0.00561002 0.0015471 0.0001115 0.0026660 0 0 0;...
                0.0102498 0.0004360 -0.0021919 -0.0010761 0.0021356 0 0;...
                -0.0013767 0.0005161 0.0005464 0.0004250 0.0013801 0.0004869 0;...
                -0.0024871 0.0011514 0.0007281 0.0013361 -0.0005090 0.0001236 0.00071695];
    
    S_Bennu = [0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            0 0 -0.0000197 0 0 0 0;...
            0 0.0015368 0.0000653 -0.0009332 0 0 0;
            0 0.0018562 0.0007749 0.0001024 0.0030684 0 0;...
            0 -0.0000786 -0.0012414 0.0002269 -0.0005038 0.0005241 0;...
            0 -0.0000084 -0.0002936 -0.0009706 -0.0006485 -0.0005993 -0.0000500];

    % order coeficient
    [M] = orderMatrix(C_Bennu, S_Bennu, Nc, Ns, normalized);
    
elseif(truth_values == "file")
    n = 16;
    R = 0.16E5;
    % count coeff
    [Nc, Ns] = countCoeff(n);

    file = table2array(readtable('Eros_CD_norm.txt'));
    [M] = orderFile(file, Nc, Ns, n, normalized);
end

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

function [Matrix] = orderMatrix(C, S, Nc, Ns, normalized)
    % Order C and S matrices in a 1x1 matrix
    Matrix = ones(Nc+Ns, 1) * NaN;

    % Point mass solution
    Matrix(1, 1) = 1;
    
    n = 2;
    m = 0;
    for j = 2:Nc
        N = n + 1;
        M = m + 1;
        if(normalized == 1)
            Norm = NormFactor(n, m);
        else
            Norm = 1;
        end

        Matrix(j) = C(N, M) * Norm;
        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n +1;
        end
    end

    n = 2;
    m = 1;
    for j = Nc + 1:Ns + Nc
        N = n + 1;
        M = m + 1;
        if(normalized == 1)
            Norm = NormFactor(n, m);
        else
            Norm = 1;
        end

        Matrix(j) = S(N, M) * Norm;
        if(m < n)
            m = m + 1;
        else
            m = 1;
            n = n + 1;
        end
    end

end

function [Matrix] = orderFile(file, Nc, Ns, degree, normalized)
    % order file given by Dan
    Matrix = ones(Nc + Ns, 1) * NaN;

    Matrix(1, 1) = file(1, 1);

    % reshape file values
    n_row = length(file(:, 1));
    n_col = length(file(1,:));

    values = reshape(file', [n_row * n_col, 1]);
    values = [1; values(5:end)];
    % index initialitation
    m = 0;
    n = 2;
    j = 2;
    kc = 2;
    ks = Nc + 1;

    while n<=degree
        if(normalized == 0)
            Norm = NormFactor(n, m);
        else
            Norm = 1;
        end

        if(m == 0)
            Matrix(kc) = values(j)*Norm;
            j = j + 1;
        else
            Matrix(kc) = values(j) * Norm;
            Matrix(ks) = values(j + 1) * Norm;
            j = j + 2;
            ks = ks + 1;
        end
        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
        kc = kc + 1;
    end
end