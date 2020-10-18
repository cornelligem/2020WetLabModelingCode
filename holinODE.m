function dxdt= holinODE(t, x, L, d)
% x(1)= antiholin
% x(2)= holin
% x(3)= holin-antiholin dimer
% L= lactate concentration
% d= degradation constant for holin
f = 23.322 / (479.3 + 7311);
basA = 479.3; %basal expression (based on fluorescence) [uM/min]
KdA = 1075; %dissociation constant [uM]
n = 1.326; %cooperativity constant
a = 7311; % alpha 1 [uM/min]
dA = 0.13;
dAT = 0.0348; % BAD, NOT GIVEN BY IGEM WIKI
kf = 0.0072; % [1/uM*min], multiplied by 60
kb = 0.0018; %[1/min], multiplied by 60
basT = (23.22 * (1303/ 2547)); %[uM/min] J23108 promoter: Multiplied by relative weakness to J23100

dxdt= zeros(3,1);

dxdt(1) = f * (basA + a.*((L / KdA )^n / (1 + (L / KdA)^n)) - (dA * x(1))) - kf * x(1) * x(2) + kb * x(3);
dxdt(2) = basT - d * x(2) - kf * x(1) * x(2) + kb * x(3);
dxdt(3) = kf * x(1) * x(2) - kb * x(3) - dAT * x(3);

end