close all;
clear all;
clc;

B = 0.1:0.05:2;
M = 10000;
h = 0.9;

EL = zeros(1,M);
sqEt = zeros(1,length(B));
Et = zeros(1,length(B));

for i = 1:length(B)
    fun = @(x) exp(-2*B(i)*x.^2);
    A = integral(fun,0,1);
    
    x = 1.0;
    count = 0;
    for j = 1:M
        if x-h/2 > 0 
            low = x - h/2;
        else
            low = 0;
        end
        if x+h/2 < 1 
            high = x + h/2;
        else
            high = 1;
        end
        xT = low + (high-low)*rand(1,1);
        
        F0 = exp(-2*B(i)*x^2)/A;
        FT = exp(-2*B(i)*xT^2)/A;
        
        r = FT/F0;
        
        if r >= 1
            x = xT;
            count = count+1;
        else
            n = rand(1,1);
            
            if n < r
                x = xT;
                count = count+1;
            end
        end
        
        EL(j) = B(i) + (0.5 - 2*B(i)^2)*x^2;
    end
    
    Et(i) = sum(EL)/M;
    sqEt(i) = sum(EL.^2)/M;
    acc = count/M
end

plot(B, Et, '.')
var = sqEt - Et.^2;

figure(2)
plot(B, var, '.')

figure(3)
y = linspace(-6,6,100);
plot(y, exp(-0.5*y.^2))
    