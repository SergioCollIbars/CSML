clear;
clc;
close all;

%%
syms p
syms x
gap = deg2rad(40);    % [rad]
n = 5;
m = 5;
phi0 = -pi/2 + gap;
phif = pi/2 + gap;

x0 = sin(phi0);
xf = sin(phif);

expression = (x*x -1)^m * diff((x*x -1)^n, n+m);
I_analytical = double(int(expression, x, x0, xf));
I_numerical = 0;
for k = 0:n
    b = sin(phif)^(2*n - 2*k +1);
    c = sin(phi0)^(2*n - 2*k +1);
    a = nchoosek(n, k) * (-1)^k / (2*n - 2*k+1) * (b - c);
    I_numerical = I_numerical + a;
end
I_numerical = factorial(2*n) * I_numerical; 

L1 = diff((x^2 - 1)^n, n) * diff((x^2 - 1)^m, m-1);
L1 = double(subs(L1, x, xf) - subs(L1, x, x0));

L2 = 0;
for k = 1:n-1
    nn = n + k;
    mm = m - 1 - k; 
    a = diff((x^2 - 1)^n, nn) * diff((x^2 - 1)^m, mm);
    L2 = L2 + (-1)^k * double(subs(a, x, xf) - subs(a, x, x0));
end


expression = (x^2 -1);
Pn = 1 / (2^n * factorial(n)) * diff(expression^n, n);
Pm = 1 / (2^m * factorial(m)) * diff(expression^m, m);

I_true = double(int(Pn*Pm, x, x0, xf));
I_computed = ((-1)^n*I_numerical + L1 + L2) * 1 / (factorial(n)^2 * 2^(2*n));
I_computed2 = (L1 + L2) * 1 / (factorial(n)*factorial(m) * 2^(n+m));

% Tested through here

%%
Nmax = 6;
H = zeros(Nmax, Nmax);

expression = (x^2 -1);
for n = 1:Nmax
    for m =1:Nmax

        Pn = 1 / (2^n * factorial(n)) * diff(expression^n, n);
        Pm = 1 / (2^m * factorial(m)) * diff(expression^m, m);

        I_true = double(int(Pn*Pm, x, x0, xf));
        H(n,m) = I_true;
    end
end
h = diag(diag(H));
H2 = H'*H;
h2 = h'*h;


figure()
subplot(1, 2, 1)
heatmap(inv(H2))

subplot(1, 2, 2)
heatmap(inv(h2))
