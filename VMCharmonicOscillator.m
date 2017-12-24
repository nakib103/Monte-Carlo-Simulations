%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREATED BY : SYED NAKIB HOSSAIN %%%
%%%%%%%% DATE : 22 DEC 2017 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% FINDING GROUND STATE FOR SIMPLE HARMONIC OSCILLATOR 
%%% USING VARIATIONAL MONTE CARLO METHOD

%%% THE RANDOM NUMBER IS GENERATED USING AN ENSAMBLE OF
%%% 50 SET WHICH CONTAIN 500 NUMBERS WITH A DISTRIBUTION
%%% THAT REFLECT THE DISTRIBUTION OF PSI(X)^2

%%% IMPORTANT LESSON TAKEN FROM THIS PROJECT:
%%% 1. THE NUMBER OF RANDOM NUMBERS MUST BE LARGE ENOUGH
%%% INITIALLY STARTED WITH 100 NIMBERS IN EACH SET BUT IT
%%% DOES NOT RESULTED IN SAME MINIMA FOR ENERGY ANS SD
%%% 2. THE RANDOM NUMBER MUST BE GENERATED WITH A 
%%% DISTRIBUTION THAT REFLECT THE IMPORTANCE FUNCTION


close all;
clear all;
clc;

B = 0:0.05:2;
M = 500;
Nstep = 50;
h = 0.9;

var = zeros(1,length(B));
Emean = zeros(1,length(B));
X = zeros(1, M);
seed = 57721566.;

for i = 1:M
    [X(i), seed] = random(0, seed);
end
mean(X)
std(X)
X = randn(1, M);
[n, xout] = hist(X, 100);
figure(3)
bar(xout, n)

for i = 1:length(B)
    EL = 0;
    ELsq = 0;
    count = 0;
    for k = 1:Nstep
        for j = 1:M
            PSIx = exp(-B(i)*X(j)^2);
            [R, seed] = random(0, seed);
            y = X(j) + (h * (R - 0.5));
            PSIy = exp(-B(i)*y^2);
        
            A = (PSIy ^ 2) / (PSIx ^ 2);
            if A > 1
                X(j) = y;
                count = count+1;
            else
                [R, seed] = random(0, seed);
                if R < A
                    X(j) = y;
                    count = count+1;
                end
            end
            
            EX = B(i) + (0.5 - 2*B(i)^2)*X(j)^2;
            EL = EL + EX;
            ELsq = ELsq + (EX * EX);
        end
    end
    
    Emean(i) = EL / (M * Nstep);
    ELsqmean = ELsq / (M * Nstep);
    var(i) = sqrt(abs((Emean(i)^ 2) - ELsqmean)) /  (M * Nstep - 1);
    acc = count/ (M * Nstep)
end

figure(1)
subplot(211)
plot(B, Emean, '.')
subplot(212)
plot(B, var, '.')

figure(2)
y = linspace(-6,6,100);
plot(y, exp(-0.5*y.^2))
    