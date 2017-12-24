close all;
clear all;
clc;

B = 0.1:0.0005:0.2;
M = 10000;
h = 0.9;

EL = zeros(1,M);
Et = zeros(1,length(B));
sqEt = zeros(1,length(B));

for i = 1:length(B)
    fun = @(r) exp(-B(i)*r.^2);
    A = integral(fun,0,1);
    
    r = 1.0;
    count = 0;
    for j = 1:M
        if r-h/2 > 0 
            low = r - h/2;
        else
            low = 0;
        end
        if r+h/2 < 1 
            high = r + h/2;
        else
            high = 1;
        end
        rT = low + (high-low)*rand(1,1);
        
        F0 = exp(-2*B(i)*r^2)/A;
        FT = exp(-2*B(i)*rT^2)/A;
        
        p = FT/F0;
        
        if p >= 1
            r = rT;
            count = count+1;
        else
            n = rand(1,1);
            
            if n < p
                r = rT;
                count = count+1;
            end
        end
    
        EL(j) = 3*B(i) - 2*B(i)^2*r^2 + 1/r;
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