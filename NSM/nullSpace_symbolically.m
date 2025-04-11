clear;
clc;
close all;

% Define symbolic variables
syms GM x y z

% Define r
r = sqrt(x^2 + y^2 + z^2);

% Construct the 5x3 matrix with symbolic entries
matrix = -[ (15*GM*x^3)/r^7 - (9*GM*x)/r^5, (15*GM*x^2*y)/r^7 - (3*GM*y)/r^5, (15*GM*x^2*z)/r^7 - (3*GM*z)/r^5;
            (15*GM*y*x^2)/r^7 - (3*GM*y)/r^5, (15*GM*x*y^2)/r^7 - (3*GM*x)/r^5, (15*GM*x*y*z)/r^7;
            (15*GM*z*x^2)/r^7 - (3*GM*z)/r^5, (15*GM*x*y*z)/r^7, (15*GM*x*z^2)/r^7 - (3*GM*x)/r^5;
            (15*GM*x*y^2)/r^7 - (3*GM*x)/r^5, (15*GM*y^3)/r^7 - (9*GM*y)/r^5, (15*GM*y^2*z)/r^7 - (3*GM*z)/r^5;
            (15*GM*x*y*z)/r^7, (15*GM*y^2*z)/r^7 - (3*GM*z)/r^5, (15*GM*z^2*y)/r^7 - (3*GM*y)/r^5];

nullSpace = null(transpose(matrix));

% Orthonormalize the basis using Gram-Schmidt
orthonormalNullSpace = orth(nullSpace);

% Display the result
disp('Orthonormal Null Space:');
disp(orthonormalNullSpace);

