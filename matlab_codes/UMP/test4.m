clear;
clc;

% Integral testing
syms x phi;
n = 2;

% 0th order LP
pol1 = (x^2 - 1)^n;
pol2 = (1- (x^2));
P_x = sqrt(2*n + 1) / (2^n * factorial(n)) * diff(pol1, n);

P_phi = subs(P_x, x, sin(phi));

I1 = int(P_x*P_x, -1, 1);
I2 = int(P_phi*P_phi*cos(phi), -pi/2, pi/2);

% 1st order LP
dP_x = sqrt(2*n + 1) / (2^n * factorial(n)) * diff(pol1, n+1);
dP_phi = diff(P_phi, phi);

I3 = int((1 - x^2) * dP_x * dP_x, -1, 1);
I4 = int(dP_phi*dP_phi*cos(phi), -pi/2, pi/2);

% 2nd order LP
ddP_x = sqrt(2*n + 1) / (2^(n) * factorial(n)) * (-x * diff(pol1, n+1) + ...
    (1-x^2)*diff(pol1, n+2));
ddP_phi = diff(dP_phi, phi);

A = sqrt(2*n + 1) / (2^n * factorial(n));
I5 = int(ddP_x*ddP_x, -1, 1);

I6 = int(ddP_phi*ddP_phi*cos(phi), -pi/2, pi/2);
I7 = sqrt(2*n + 1) / (2^n * factorial(n)) *  sqrt(2*n + 1) / (2^n * factorial(n))* (x^2 * diff(pol1, n+1)*diff(pol1,n+1) - ...
    n*(n+1)*(1-x^2)*diff(pol1, n+2)*diff(pol1,n));
I7 = int(I7, -1, 1);
I8 = int(ddP_x * P_x,-1,1);
ex = x*diff(pol1,n+1)*diff(pol1, n) - n*(n+1)*diff(pol1, n)^2;
I9 = int(A*A*ex,-1, 1);
d = x*diff(pol1,n+1)*diff(pol1, n);
d = int(d, -1, 1);


A = 1 / (factorial(n)^2 * 2^(2*n));
n2 = factorial(2*n);
nf = factorial(n);

I2 = int(pol2* diff(pol1, n)*diff(pol1, n+2),-1,1);
I2_sol = -n2*nf^3/factorial(n-2)*2^(1+2*n)/factorial(2*n +1);
uv = x^2*diff(pol1,n+1)*diff(pol1, n);
uv  = subs(uv, 1) - subs(uv, -1);
uv_sol = nf* 2^(2*n)*n*factorial(n+1);
I11 = diff(pol1, n) * diff(x^2*diff(pol1, n+1));
I11 = int(-I11, -1, 1);
I11_sol = -n2/factorial(n-1)*factorial(n+1)*nf^2*2^(1+2*n)/factorial(2*n+1);

I1  = int(x^2*diff(pol1,n+1)*diff(pol1, n+1),-1,1);
I1_sol = uv_sol + I11_sol;
a = sqrt(2*n + 1) / (2^(n) * factorial(n));
I = a*a * (-n*(n+1)*I2_sol + I1_sol);

% 2nd order LP split integrals

B = (-1)^n*factorial(n)^2*(2^(1+2*n))/factorial(2*n+1);
C = (factorial(n)^2*2^(2*n));
sol2 = factorial(n)*factorial(2*n)*B/factorial(n-2);
exp2 = (1-x^2)*diff(pol1,n)*diff(pol1,n+2);
exp1 = (x^2)*diff(pol1, n+1)*diff(pol1, n+1);
sol1 = factorial(n+1)*factorial(2*n)*B/factorial(n-1);
d = factorial(n)*n*factorial(n+1)*2^(2*n);
e = factorial(n+1)*factorial(2*n)/factorial(n-1)*B + d;

a = diff(pol1, n+1);
b = diff(pol1, n);
c = subs(a*b, 1) - subs(a*b,-1);


i2 = int(exp2, -1, 1);
i1 = int(exp1, -1, 1);
r = (x^2) * diff(pol1, n+1);
u = diff(r, n+1);
A = factorial(2*n) * factorial(n+1) / factorial(n-1);

I = (2*n +1) * B / C * factorial(2*n) * (factorial(n) / factorial(n-1) - n*(n+1));





