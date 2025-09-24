clear;
clc;
close all;
format long g;
addpath("data/")
addpath("functions/")

%%          JPL WEBSITE CR3BP ORBIT FAMILES CONVERTER
% Date: 09/16/2025
% Description: Parser the information from the .csv files downloaded from
% the JPL website and store them in a list formating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_name  = "periodic_orbits";          % . csv file
output_name = "JPL_EM_Vert_L2_Family";   % . mat file

% Load CSV into a table
T = readtable(input_name + ".csv"); 

N = height(T);                 % number of rows
trajFam = cell(N,1);           % N×1 cell array

disp('Reading and parsing file ...')
for i = 1:N
    % Create a struct with fields equal to the column names
    rowStruct = struct();
    rowStruct.("mu") = 0.0121505853505625;
    rowStruct.("iState") = T{i, 2:7}';
    rowStruct.("period") = T{i, 9};
    rowStruct.("CJ") = T{i, 8};

    trajFam{i} = rowStruct;   % store in the cell
end
disp('      DONE!');
save(output_name+".mat", "trajFam");