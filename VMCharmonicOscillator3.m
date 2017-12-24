%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREATED BY : SYED NAKIB HOSSAIN %%%
%%%%%%%% DATE : 22 DEC 2017 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% FINDING GROUND STATE FOR SIMPLE HARMONIC OSCILLATOR 
%%% USING VARIATIONAL MONTE CARLO METHOD

%%% WE GENERATE A SET OF 50 RANDOM NUMBERS EACH CONTAINING 
%%% 10000 RANDOM NUMBERS USING METROPLOIS ALGORITHM 

%%% WHY DOES THE ENERGY AT HIGH LAMDA IS LOWER THAN AT 0.5?


close all;
clear all;
clc;

B = 0.1:.05:2;
M = 10000;
Nstep = 50;
h = 0.9;

var = zeros(1,length(B));
Emean = zeros(1,length(B));

for i = 1:length(B)
    EL = 0;
    ELsq = 0;
    count = 0;
    for k = 1:Nstep
        x = 0;
        for j = 1:M
            low = x - h / 2;
            high = x + h / 2;
            xprime = low + (high - low) * rand();
            
            F = exp(- B(i) * x ^ 2);
            Fprime = exp(- B(i) * xprime ^ 2);
            
            A = Fprime / F;
            if A > 1
                x = xprime;
            elseif A > rand()
                x = xprime;
            end
            
            EX = B(i) + (0.5 - 2*B(i)^2)*x^2;
            EL = EL + EX;
            ELsq = ELsq + EX ^ 2;
        end        
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
    