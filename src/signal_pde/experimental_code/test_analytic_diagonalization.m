clear
clc
close all hidden

syms D xsi k_on k_off lambda A t real

A = [-D*xsi^2-k_on, k_off ; k_on, -k_off]*t;

[PP, DD] = eig(A); 