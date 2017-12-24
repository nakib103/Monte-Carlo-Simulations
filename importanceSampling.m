%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREATED BY : SYED NAKIB HOSSAIN %%%
%%%%%%%% DATE : 23 DEC 2017 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% IMPLEMENTING IMPORTANCE SAMPLING IN MONTE CARLO
%%% IN INTEGRATION. THE FUNCTION CONSIDERED ARE
%%% 1 / (X^2 + COS(X)^2) AND EXP(-X) / (1 + (X - 1)^2)

%%% THE IMPORTANCE FUNCTION FOR THE FUNCTIONS ARE 
%%% CONSIDERED TO BE, A*EXP(- LAMDA * X) IN BOTH CASES.
%%% LAMDA IS OPTIMIZED TO APPROACH THE INTEGRAL VALUE.

%%% THE FIRST FUNCTION IS ALSO DONE IN CRUDE MC METHOD

close all;
clear all;
clc;

%% CRUDE MC FOR FIRST FUNCTION

N = 10000;
Nstep = 50;

Favg = zeros(1, Nstep);
for i = 1:Nstep
    X = pi*rand(1, N);
    
    F = pi*sum(1.0 ./ (X.^2 + cos(X).^2));
    Favg(i) = F / N;
    
    Fsq = sum((pi*(1.0 ./ (X.^2 + cos(X).^2))) .^ 2) / N;
    Var = sqrt(Fsq - Favg);
end

mean(Favg)
std(Favg)

%% IMPORTANCE SAMPLING MC FOR FIRST FUNCTION

lamda = 0.05:0.05:1.6;
Imean = zeros(Nstep, length(lamda));
Var = zeros(Nstep, length(lamda));
X = rand(Nstep, N);

for k = 1:length(lamda)
    
    r = -(1 / lamda(k)) * log(1 - ((1 - exp(-lamda(k)*pi)) * X));
    %r = -(1 / lamda(k)) * log(((1 - exp(-lamda(k)*X)) / lamda(k)));
    for i = 1:Nstep
        I = 0;
        Isq = 0;
        for j = 1:N
            G = (lamda(k) / (1 - exp(-lamda(k)*pi)) ) * exp(-lamda(k)*r(j));
            F = (1.0 ./ (r(j)^2 + cos(r(j))^2));
            EX = F / G;

            I = I + EX;
            Isq = Isq + (EX * EX);
        end
    
    Imean(i, k) = I / N;

    Isqmean = Isq / N;
    Var(i, k) = abs(Isqmean - Imean(k)^2);
    end
end

mean(Imean(find(lamda == 0.8), :))
std(Imean(find(lamda == 0.8), :))

subplot(211)
plot(lamda, mean(Var))
subplot(212)
plot(lamda, mean(Imean))

%% IMPORTANCE SAMPLING MC FOR SECOND FUNCTION

lamda = 0.05:0.05:2.5;
Imean = zeros(Nstep, length(lamda));
Var = zeros(Nstep, length(lamda));
X = rand(Nstep, N);

for k = 1:length(lamda)
    r = - (1 / lamda(k)) * log(1 - X);
    for i = 1:Nstep
        I = 0;
        Isq = 0;
        for j = 1:N
            F = exp(-r(j)) / (1 + (r(j) - 1) ^ 2);
            G = lamda(k) * exp(- lamda(k) * r(j));
            EX = F / G;
            
            I = I + EX;
            Isq = Isq + (EX * EX);
        end
        
        Imean(i, k) = I / N;
        
        Isqmean = Isq / N;
        Var(i, k) = abs(Isqmean - Imean(k)^ 2);
    end
end

figure(2)
subplot(211)
plot(lamda, mean(Var))
subplot(212)
plot(lamda, mean(Imean))