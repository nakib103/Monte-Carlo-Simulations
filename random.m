function [r, seed] = random(X, seed)
    A = 3141592621.;
    C = 10000000000.;
    P = 2718281829.;
    
    R = A * seed + C;
    seed = mod(R, P);
    r = seed / P + X;
end