clear
clc
close all hidden

% syms D xsi k_on k_off lambda A t real
% 
% A = [-D*xsi^2-k_on, k_off ; k_on, -k_off]*t;
% 
% [PP, DD] = eig(A);


D = 200 + 500 * rand();
xsi = 256 * rand();
k_on = 0.5 + 3 *rand();
k_off = 1 + 5 *rand();
t = 25 * rand();

A = [-D*xsi^2-k_on, k_off ; k_on, -k_off]*t
expmA = expm(A)

PP = [ -(k_on - k_off + (D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2 + k_on^2 + 2*k_on*k_off + k_off^2)^(1/2) + D*xsi^2)/(2*k_on), -(k_on - k_off - (D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2 + k_on^2 + 2*k_on*k_off + k_off^2)^(1/2) + D*xsi^2)/(2*k_on) ; ...
                                                                                                                                 1,                                                                                                                             1];
 
DD = [ -(t*(k_on + k_off + (D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2 + k_on^2 + 2*k_on*k_off + k_off^2)^(1/2) + D*xsi^2))/2,                                                                                                                        0 ; ...
                                                                                                                              0, -(t*(k_on + k_off - (D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2 + k_on^2 + 2*k_on*k_off + k_off^2)^(1/2) + D*xsi^2))/2  ];

PPinv = [ -k_on/(2*k_on*k_off + k_on^2 + k_off^2 + D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2)^(1/2), -(k_on - k_off - (D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2 + k_on^2 + 2*k_on*k_off + k_off^2)^(1/2) + D*xsi^2)/(2*(2*k_on*k_off + k_on^2 + k_off^2 + D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2)^(1/2)) ; ...
         k_on/(2*k_on*k_off + k_on^2 + k_off^2 + D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2)^(1/2),  (k_on - k_off + (D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2 + k_on^2 + 2*k_on*k_off + k_off^2)^(1/2) + D*xsi^2)/(2*(2*k_on*k_off + k_on^2 + k_off^2 + D^2*xsi^4 + 2*D*k_on*xsi^2 - 2*D*k_off*xsi^2)^(1/2))];
 

% PP
% DD
% PPinv

Aprim = PP*DD*PPinv
expmAprim = PP*[exp(DD(1,1)) 0 ; 0 exp(DD(2,2))]*PPinv


sum( (expmA(:) - expmAprim(:)).^2 )