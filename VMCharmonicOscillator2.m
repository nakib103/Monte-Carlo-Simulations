%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREATED BY : SYED NAKIB HOSSAIN %%%
%%%%%%%% DATE : 22 DEC 2017 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% FINDING GROUND STATE FOR SIMPLE HARMONIC OSCILLATOR 
%%% USING VARIATIONAL MONTE CARLO METHOD

%%% THE RANDOM NUMBER IS GENERATED USING AN ENSAMBLE OF 50 
%%% SET WHICH CONTAIN 500 NUMBERS WITH A DISTRIBUTION
%%% THAT REFLECT THE DISTRIBUTION OF PSI(X)^2. WE GENERATE
%%% RANDOM NUMBERS USING DEFAULT NORMRND() FUNCTION IN MATLAB

%%% WHY DOES THE ENERGY CCURVE IS CONVEX ISTEAD OF CONCAVE? 


close all;
clear all;
clc;

rng('shuffle')

B = 0:0.05:2;
M = 500;
Nstep = 50;
h = 0.9;

var = zeros(1,length(B));
Emean = zeros(1,length(B));

X = normrnd(0, 1, 1, M);
mean(X)
std(X)
[n, xout] = hist(X, 100);
figure(1)
%bar(xout, n)

for i = 1:length(B)
    EL = 0;
    ELsq = 0;
    count = 0;
    for k = 1:Nstep
        for j = 1:M
            PSIx = exp(-B(i)*X(j)^2);
            R = normrnd(0, 1);
            Y = X(j) + (h * (R - 0.5));
            PSIy = exp(-B(i)*Y^2);

            A = (PSIy ^ 2) / (PSIx ^ 2);
            if A > 1
                X(j) = Y;
            else
                r = normrnd(0, 1);
                if A > r
                    X(j) = Y;
                end
            end
        end

        EX = B(i) + (0.5 - 2*B(i)^2)*X(j).^2;
        EL = sum(EX);
        ELsq = sum(EX.^2);
    end
    
    
    Emean(i) = EL / (M * Nstep);
    ELsqmean = ELsq / (M * Nstep);
    var(i) = sqrt(abs((Emean(i)^ 2) - ELsqmean)) /  (M * Nstep - 1);
    %acc = count/ (M * Nstep)
end

figure(2)
subplot(211)
plot(B, Emean, '.')
subplot(212)
plot(B, var, '.')

figure(3)
y = linspace(-6,6,100);
plot(y, exp(-0.5*y.^2))
    