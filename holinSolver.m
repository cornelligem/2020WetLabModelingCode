deg= 0.0348; % DEGRADATION IS SOMETHING WE COULD MODIFY
L = [0.01, 0.1, 1, 10, 100] * 1000;
options = odeset('RelTol', 1e-10);
tspan = [0 300];
y0 = [ 0; 0; 0 ];
fprintf("Steady state toxin level is: \n");

figure;

for i= 1:length(L)
[t, x] = ode45(@(t,x) holinODE(t,x,L(i),deg), tspan, y0);
subplot(1, length(L), i)
plot(t, x, 'LineWidth', 1.5)
legend(["antiholin", "holin", "dimer"]);
title("Lactate = " + L(i) / 1000 + " mM");
ylabel("Species concentration (uM)")
xlabel("Time (min)")
ylim([0 500])
fprintf("[Lactate] = %f mM: %f molecules\n", L(i) / 1000, getMolecule(x(end, 2)));
end
%%

D= [0.0001, 0.01, 1, 100];
for i= 1:length(D)
[t1, x1] = ode45(@(t,x) holinODE(t,x,1 * 1000,D(i)), tspan, y0, options);
[t10, x10] = ode45(@(t,x) holinODE(t,x,10 * 1000,D(i)), tspan, y0, options);
subplot(1, length(D), i)
hold on
plot(t1, x1(:,2), '-or', 'LineWidth', 1.5)
plot(t10, x10(:,2), '-*b', 'LineWidth', 1.5)
hold off
legend(["toxin with L= 1", "toxin with L= 10"]);
title("Degradation rate = " + D(i) + " min^-^1");
ylabel("Species concentration (uM)")
xlabel("Time (min)")
fold_increase= x1(end, 2) / x10(end, 2);
fprintf("Degradation rate = %f per min: %f fold decrease\n", D(i), fold_increase);
fprintf(2, "    [L] = 1: %f  |  [L] = 10: %f\n", getMolecule(x1(end, 2)), getMolecule(x10(end, 2)));
end

%% Will take a while (5-10 minutes), be patient
% Analyze 3D
L_vec= linspace(0, 30, 15);
D_vec= logspace(-4, 4, 15);
output= zeros(length(L_vec), length(D_vec));
for i = 1:length(L_vec)
    for j= 1:length(D_vec)
        [t, x] = ode45(@(t,x) holinODE(t,x,L_vec(i) * 1000, D_vec(j)), tspan, y0, options);
        output(i,j) = getMolecule(x(end, 2));
        disp(i + ", " + j);
    end
end
%%
figure;
imagesc(L_vec, D_vec, output);
xlabel("Lactate Concentration (mM)");
ylabel("Degradation constant (min^-^1)");
colorbar;
caxis([0, 250])

%% 

function num_mol= getMolecule(uM)
    vol_E_coli= 0.6e-15;
    mole= 6.022e23;
    num_mol = (uM / 1e6) * (vol_E_coli * mole);
end
