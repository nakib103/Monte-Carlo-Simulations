%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREATED BY : SYED NAKIB HOSSAIN %%%
%%%%%%%% DATE : 23 DEC 2017 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% IMPLEMENTING METROPOLIS ALGORITHM IN SAMPLING RANDOM
%%% NUMBER WITH DISTRIBUTION LIKE THE IMPORTANCE FUNCTION

%%% WE FIRST SAMPLE FROM THE FUNCTION EXP(- 0.5 * X ^ 2)
%%% FOR PRACTICE. THEN SAMPLE THE FUNCTION A * EXP(X ^ 2 - 1)  
%%% AND USE IT TO INTEGRATE THE FUNCTION X * SQRT(1 - X ^ 2)

close all;
clear all;
clc;

%% PRACTICE: SAMPLING USING METROPOLIS ALGORITHM 
M = 100000;
h = 5;
x = 0;
X = zeros(1, M);

count = 0;
for i = 1:M
    low = x - h / 2;
    high = x + h / 2;
    xprime = low + (high - low)*rand();
    
    F = exp(- 0.5 * x ^ 2);
    Fprime = exp(- 0.5 * xprime ^ 2);
    
    A = Fprime / F;
    if A > 1
        x = xprime;
        count = count + 1;
    elseif A > rand()
        x = xprime;
        count = count + 1;
    end
    
    X(i) = x;
end 
disp(count / M)
[y, n] = hist(X, 50);
figure(1)
bar(n, y);

%% INTEGRATION USING IMPORTANCE SAMPLING MONTE CARLO AND METROPOLIS ALGORITHM

M = 100000;
Nstep = 50;
h = 0.42;
x = 0;
X = zeros(Nstep, M);
I = zeros(1, Nstep);

for j = 1:Nstep
    count = 0;
    for i = 1:M
        low = max(0, x - h / 2);
        high = min(1, x + h / 2);
        xprime = low + (high - low)*rand();

        F = 2.1615*( exp(x ^ 2) - 1);
        Fprime = 2.1615*( exp(xprime ^ 2) - 1);

        A = Fprime / F;
        if A > 1
            x = xprime;
            count = count + 1;
        elseif A > rand()
            x = xprime;
            count = count + 1;
        end

        X(j, i) = x;
    end 
    F = X(j, :) .* sqrt(1 - X(j, :) .^ 2);
    G = 2.1615 .* (exp(X(j, :) .^ 2) - 1);
    I(j) = sum(F ./ G) / M;
end

disp(count / M)
[y, n] = hist(X(1, :), 50);
figure(2)
bar(n, y);

mean(I)
std(I)