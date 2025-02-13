clear all;
clc;
close all
%%          LEGENDRE POLINOMIALS TEST

syms lambda phi 
n_max = 10;
normalized = 1;

grid_ddP_ddP = meshgrid(1:n_max, 1:n_max)*0;
grid_ddP_P = meshgrid(1:n_max, 1:n_max)*0;

Nc = 0;
for j = 1:n_max+1
    Nc  = Nc + j;
end
Nc = Nc - 2;
Ns = Nc - n_max;
Nx = Nc + Ns;

% define out matrices
output = zeros(n_max^2, 5);
p = 1;
for j = 1:n_max % n1
    n1 = j;
    m1 = 0;
    for k = 1:n_max
        clc;
        disp('n1 = ' + string(j) + ' n2 = ' + string(k))
        n2 = k;
        m2 = 0;

        P2 = assocLegendre(n2, m2);
        dP1 = Diff_assocLegendre(n1, m1);
        dP2 = Diff_assocLegendre(n2, m2);
    
        ddP1 = Diff2_assocLegendre(n1, m1);
        ddP2 = Diff2_assocLegendre(n2, m2);
    
        if(normalized)
            P2 = NormFactor(n2, m2) * P2;
            dP1 = NormFactor(n1, m1) * dP1;
            dP2 = NormFactor(n2, m2) * dP2; 
            ddP1 = NormFactor(n1, m1) * ddP1;
            ddP2 = NormFactor(n2, m2) * ddP2;
        end
        Y2 = P2 * cos(m2 * lambda);
        dY1 = dP1 * cos(m1 * lambda);
        dY2 = dP2 * cos(m2 * lambda);
    
        ddY1 = ddP1 * cos(m1 * lambda);
        ddY2 = ddP2 * cos(m2 * lambda);
    
        expression = dY1 * dY2 * cos(phi);
    
        expression = int(expression, phi, -pi/2, pi/2);
        solution = int(expression, lambda, 0, 2*pi);
        output(p, 3) = double(solution);
        
        expression = ddY1 * ddY2 * cos(phi);
        
        solution = int(expression, phi, -pi/2, pi/2);
        output(p, 4) = double(solution);
        
        % special expression
        expression = ddY1 * Y2 * cos(phi);
        
        solution = int(expression, phi, -pi/2, pi/2);
        output(p, 5) = double(solution);

        output(p, 1) = n1;
        output(p, 2) = n2;

        %fill grid
        if(k <= j)
            grid_ddP_ddP(j, k) = output(p, 4);
            grid_ddP_P(j, k) = output(p, 5);
        end

        % update counter
        p = p + 1;
    end
end

% create output table
varNames = ["n1","n2","dY_n1 * dY_n2", "ddY_n1 * ddY_n2", "special expression"];
tab = array2table(output(:, 1:5),'VariableNames',varNames);
disp(tab);

%% PLOTS

coeffVal = {};
for j = 1:n_max
    coeffVal{end+1} = string(j);
end

% plot tendency
xAxis = 1:n_max;
data = zeros(length(xAxis), 2);
n = 1;
for j = 1:length(output(:, 1))
    if(output(j, 3) ~= 0)
        data(n, 1) = output(j, 3);
        data(n, 2) = output(j, 4);
         n = n+1;
    end
end
figure()

trend = (xAxis + 1).*xAxis;
plot(xAxis, data(:, 1)./(4*pi), '*')
hold on;
plot(xAxis, trend, '--');
xlabel('harmonic degree, n');
ylabel('value / 4\pi')
title('integral first derivative LP product: |dY_{n0}|^2');

figure()

trend = (xAxis + 1).*xAxis.*(xAxis.^2 - 0.5);
trend(1) = 0;
plot(xAxis, data(:, 2)./(4*pi), '*')
hold on;
plot(xAxis, trend, '--');
xlabel('harmonic degree, n');
ylabel('value / 4\pi')
title('integral second derivative LP product: |ddY_{n0}|^2');

figure();
plot(linspace(1, n_max, n_max^2), output(:, 4)./(4*pi), '*')
hold on;
plot(xAxis, trend, '--');
title('|ddY_{n0}|^2')

figure();
trend = xAxis.^2;
trend(1) = 0;
plot(linspace(1, n_max, n_max^2), -output(:, 5)./(4*pi), '*')
hold on;
plot(xAxis, trend, '--');
title('ddY_{n0} * Y_{n0}')

% plot surface
figure();
X = linspace(1, n_max, n_max);
Y = linspace(1, n_max, n_max);
heatmap(X, Y, grid_ddP_ddP)
title('Integral value of the product between 2nd order LP')

figure();
X = linspace(1, n_max, n_max);
Y = linspace(1, n_max, n_max);
heatmap(X, Y, grid_ddP_P)
title('Integral value of the product between 2nd and 0th order LP')


%%              FUNCTIONS
function plm = assocLegendre(l,m)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l,x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm = ((1 - x^2)^(m/2))*legDiff;

    plm = subs(plm, x, sin(phi));
end

function dPlm = Diff_assocLegendre(l,m)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l,x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm = ((1 - x^2)^(m/2))*legDiff;

    plm = subs(plm, x, sin(phi));
    dPlm = diff(plm, phi);
end

function ddPlm = Diff2_assocLegendre(l,m)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l,x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm = ((1 - x^2)^(m/2))*legDiff;

    plm = subs(plm, x, sin(phi));
    dPlm = diff(plm, phi);
    ddPlm = diff(dPlm, phi);
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