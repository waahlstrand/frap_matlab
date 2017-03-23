clear
clc
close all hidden

number_of_points = 256;
[X1, X2] = meshgrid(1:number_of_points, 1:number_of_points);

X1 = X1 - mean(X1(:));
X2 = X2 - mean(X2(:));

alpha = 0.2;
f = exp( - alpha^2 * (X1.^2 + X2.^2) );

figure, imagesc(f)

F = fft2(f);
F = fftshift(F);

figure, imagesc(abs(F))
