clear;
clc;
close all;

%%                  FREQUENCY PER DEGREE ANALYTICAL
% Description: Compute the excited frequency for a certain degree in the
% gradiometer signal. Assuming a polar, circular orbit and non-rotating
% planet.
%
% Date: 10/02/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference radius & planet mass
Ref = 1737.4E3;                     % [m] 
GM  = 4.902800118E12;               % [m^3/s^2] 

% Orbit altitude
N = 5;
hmin = 1E3;                         % [m]
hmax = 3000E3;                       % [m]

separation = round((hmax - hmin) / (N-1));

H = hmin:separation:hmax;

% degree expansion
n_range = 2:1:100;

f_n = ones(N, length(n_range)) * NaN;
for k = 1:N
    r = Ref + H(k);
    f_n(k, :) = (sqrt(GM /r)/(2 * pi * Ref)).* n_range;
end

figure()
semilogy(n_range, f_n, 'LineWidth', 2)

%%
r = Ref + 3000E3;
f_orb   = sqrt(GM/r^3)/ (2*pi);
f_spin = 1/(27.3 * 86400).*1; 

f2 = grad_lines(2, 0, f_orb, f_spin);
%%
nmax   = 10;
i_deg  = 90;                % polar
f_orb  = 1/6000;            % ~100 min orbit
f_spin = 1/86164;           % 1 sidereal day

[atoms, groups] = grad_lines_map(nmax, i_deg, f_orb, f_spin);

% Inspect exactly which harmonic creates a line:
k = 12;   % pick any atom
fprintf('f=% .6e Hz  comes from  %s\n', atoms(k).f, atoms(k).label);

% Plot unique sticks with all contributors shown in the title/tooltip
hold on; xline([groups.f], '--'); hold off;

% Example: print contributors for the 5th grouped line
g = 5;
fprintf('f=% .6e Hz has contributors:\n', groups(g).f);
disp({atoms(groups(g).contributors).label}.');



%% FUNCTIONS
function f = grad_lines(n, i_deg, f_orb, f_spin)
% n: degree
% i_deg: inclination in degrees
% f_orb: orbital frequency [Hz]
% f_spin: body spin frequency [Hz]
ci = cosd(i_deg);
f_rel = abs(f_orb*ci - f_spin);

f = [];                     % list of expected line frequencies
for m = 0:n
    qset = m:2:n;           % parity selection
    for q = qset
        f = [f, abs(q*f_orb + m*f_rel), abs(q*f_orb - m*f_rel)];
    end
end
f = unique(f);              % collapse duplicates
end


function [atoms, groups] = grad_lines_map(nmax, i_deg, f_orb, f_spin, varargin)
%GRAD_LINES_MAP  Enumerate and label gravity(-gradient) spectral lines.
%
% [atoms, groups] = grad_lines_map(nmax, i_deg, f_orb, f_spin, 'Name',Value,...)
%
% Inputs
%   nmax    : maximum spherical-harmonic degree to enumerate (>=0)
%   i_deg   : orbital inclination in degrees (0=equatorial, 90=polar)
%   f_orb   : orbital frequency [Hz] = n0/(2*pi)
%   f_spin  : body spin frequency [Hz] = omega_b/(2*pi)
%
% Name-Value options
%   'Tol'         : frequency grouping tolerance (default 1e-10 Hz)
%   'IncludeDC'   : true to keep f=0 lines (default false)
%   'OnlyPlus'    : true to only make q*f_orb + m*f_rel (default false: makes both + and -)
%
% Outputs
%   atoms  : struct array with one entry per (n,m,q,combination) *before* grouping
%            fields: n,m,q,comb ('+','-','0'), f [Hz], label, latex
%   groups : struct array of unique frequencies (within Tol), with contributors[]
%
% Notes
% - Allowed latitudinal indices follow the Legendre parity rule:
%     q = n, n-2, n-4, ... down to >= m
%   (i.e., q has the same parity as n and q >= m).
% - Combination frequencies with rotation:
%     f = | q*f_orb ± m*f_rel |,  where  f_rel = | f_orb*cos(i) - f_spin |
% - For non-rotating body, use f_spin = 0.
%
% Example
%   [atoms, groups] = grad_lines_map(6, 90, 1/6000, 1/86164);
%   % Print first few atoms:
%   for k=1:10
%     fprintf('f=% .6e Hz  <- %s\n', atoms(k).f, atoms(k).label);
%   end
%   % Overlay sticks:
%   hold on; xline([groups.f], '--'); hold off;

% ---- options
p = inputParser;
p.addParameter('Tol', 1e-10);
p.addParameter('IncludeDC', false);
p.addParameter('OnlyPlus', false);
p.parse(varargin{:});
Tol       = p.Results.Tol;
IncludeDC = p.Results.IncludeDC;
OnlyPlus  = p.Results.OnlyPlus;

ci    = cosd(i_deg);
f_rel = abs(f_orb*ci - f_spin);

atoms = struct('n',{},'m',{},'q',{},'comb',{},'f',{},'label',{},'latex',{});
idx = 0;

for n = 0:nmax
    for m = 0:n
        % Latitudinal harmonics respecting parity of n: q = n, n-2, ..., >= m
        qvals = n:-2:0;
        qvals = qvals(qvals >= m);
        for q = qvals
            if m == 0
                % No longitudinal modulation term when m=0
                f0 = abs(q*f_orb);
                if IncludeDC || f0 > 0
                    idx = idx + 1;
                    atoms(idx) = make_atom(n,m,q,'0',f0);
                end
            else
                % With rotation, produce ± combinations; without rotation f_rel may be 0
                fplus  = abs(q*f_orb + m*f_rel);
                fminus = abs(q*f_orb - m*f_rel);

                if OnlyPlus
                    if IncludeDC || fplus > 0
                        idx = idx + 1; atoms(idx) = make_atom(n,m,q,'+',fplus);
                    end
                else
                    if IncludeDC || fplus > 0
                        idx = idx + 1; atoms(idx) = make_atom(n,m,q,'+',fplus);
                    end
                    if IncludeDC || fminus > 0
                        idx = idx + 1; atoms(idx) = make_atom(n,m,q,'-',fminus);
                    end
                end
            end
        end
    end
end

% Sort atoms by frequency
[~,ord] = sort([atoms.f]);
atoms = atoms(ord);

% Group near-duplicate frequencies and list contributors
groups = struct('f',{},'contributors',{});
if isempty(atoms), return; end

f_list = [atoms.f];
used = false(size(f_list));
g = 0;

for k = 1:numel(atoms)
    if used(k), continue; end
    g = g + 1;
    fk = f_list(k);
    close_idx = find(abs(f_list - fk) <= Tol);
    used(close_idx) = true;
    groups(g).f = mean(f_list(close_idx));   % representative frequency
    groups(g).contributors = close_idx(:)';  % indices into atoms
end

% ---- nested helper
    function a = make_atom(nn, mm, qq, comb, ff)
        a.n = nn; a.m = mm; a.q = qq; a.comb = comb; a.f = ff;
        a.label = sprintf('n=%d, m=%d, q=%d, %s', nn, mm, qq, combstr(comb));
        a.latex = sprintf('$n=%d,\\;m=%d,\\;q=%d,\\;%s$', nn, mm, qq, combsym(comb));
    end
    function s = combstr(c)
        switch c, case '+', s = 'q f_{orb} + m f_{rel}';
                    case '-', s = 'q f_{orb} - m f_{rel}';
                    otherwise, s = 'q f_{orb}'; end
    end
    function s = combsym(c)
        switch c, case '+', s = '+';
                    case '-', s = '-';
                    otherwise, s = '0'; end
    end
end
