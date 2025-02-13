clear;
clc;
close all;

%%      PLOT UPM RMS ERROR
% Description: Plot the total error (numerical - analytical) for every RMS
% value contained in the data files

% inputs
addpath("//Users/sergiocollibars/Desktop/uncertainty model/UPM_Bennu/")
addpath("//Users/sergiocollibars/Desktop/uncertainty model/UPM_Eros/")

% asteroid selection
asteroid = "Bennu";         % options: Bennu / Eros
% load data
if(asteroid == "Bennu")
    a = load("a_Bennu.mat").variable(1, :);
    RSS = load("RSS_Bennu.mat").MC_result;
    
    % simulation variables
    GM =  5.1355;
    W = 6.32E-13^2;
    L = 77760;
    R = 244.9400088;
elseif(asteroid == "Eros")
    a = load("a_Eros.mat").variable(1, :);
    RSS = load("sigma_Eros.mat").MC_result;
    
    % simulation variables
    GM =  4.5960443E5;
    W = 6.32E-13^2;
    L = 91584;
    R = 16000;
end
N_points = length(a(1, :));
n_max = length(RSS(:, 1));

% error variables
err = zeros(n_max, N_points);

for n = 1:n_max
    for k = 1:N_points
        % radius
        r = a(k);
        j = n;
        if(n == 1)
            n_val = 0;
        else
            n_val = n;
        end
        sigma  = compute_UMP(GM, R, r, W, L, n_val); 

        % error 
        err(j, k) = abs(RSS(j, k) - sigma);
    end
end

% plot mear error
[p,~]=numSubplots(n_max);
row = p(1);
col = p(2);

[a, I] = sort(a);

figure();
for j =1:n_max
    subplot(row, col, j);
    plot(a./1000, err(j, I), 'LineWidth', 1.5, 'Color','g')
    if(j == 1)
        val = 0;
    else
        val = j;
    end
    title('n = ' + string(val))
    xlabel('a [Km]')
    ylabel('\epsilon [-]')
end
sgtitle("RMS absolute error")





%% FUNCTION
function [sigma] = compute_UMP(GM, R, r, W, L, n)
    
    % polynomial
    pol= (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
    
    % RSS sigma value ( multiplied by scaling factor)
    sigma = ((r.^3)./GM).* ((r./R).^n) * 1/sqrt(L) *...
            sqrt(W)/sqrt(pol) * sqrt(2*n + 1);
end

function [p,n]=numSubplots(n)
    % function [p,n]=numSubplots(n)
    %
    % Purpose
    % Calculate how many rows and columns of sub-plots are needed to
    % neatly display n subplots. 
    %
    % Inputs
    % n - the desired number of subplots.     
    %  
    % Outputs
    % p - a vector length 2 defining the number of rows and number of
    %     columns required to show n plots.     
    % [ n - the current number of subplots. This output is used only by
    %       this function for a recursive call.]
    %
    %
    %
    % Example: neatly lay out 13 sub-plots
    % >> p=numSubplots(13)
    % p = 
    %     3   5
    % for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
    %
    %
    % Rob Campbell - January 2010
       
        
    while isprime(n) & n>4, 
        n=n+1;
    end
    p=factor(n);
    if length(p)==1
        p=[1,p];
        return
    end
    while length(p)>2
        if length(p)>=4
            p(1)=p(1)*p(end-1);
            p(2)=p(2)*p(end);
            p(end-1:end)=[];
        else
            p(1)=p(1)*p(2);
            p(2)=[];
        end    
        p=sort(p);
    end
    %Reformat if the column/row ratio is too large: we want a roughly
    %square design 
    while p(2)/p(1)>2.5
        N=n+1;
        [p,n]=numSubplots(N); %Recursive!
    end
end