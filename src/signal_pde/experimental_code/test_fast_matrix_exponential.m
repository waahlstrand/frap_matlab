clear
clc
close all hidden

D = 500;
k_on = 1;
k_off = 3;
A = [-D*rand()-k_on, k_off ; k_on, -k_off];

tic
for i = 1:1000
    B_MATLAB = expm(A);
end
toc

tic
for i = 1:1000
    B = [1 0;0 1];
    T = [1 0;0 1];
    for j = 1:100
        T = T * A / j;
        B = B + T;
    end
end
toc
