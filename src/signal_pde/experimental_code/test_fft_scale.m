clear
clc
close all hidden

number_of_points = 256;
[X1, X2] = meshgrid(1:number_of_points, 1:number_of_points);

X1 = X1 - mean(X1(:));
X2 = X2 - mean(X2(:));

alpha = 0.05;
f = exp( - alpha * (X1.^2 + X2.^2) );

F = fft2(f);
F = fftshift(F);

figure, imagesc(f)
figure, imagesc(abs(F))

[K1, K2] = meshgrid(1:number_of_points, 1:number_of_points);
K1 = K1 - mean(K1(:));
K2 = K2 - mean(K2(:));

K1 = K1 * 2*pi/number_of_points;
K2 = K2 * 2*pi/number_of_points;

F_analytical = sqrt(pi/alpha) * exp( - pi^2 * (K1.^2 + K2.^2) / alpha );
figure, imagesc(F_analytical)