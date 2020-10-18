%% Graph with Holin, Antiholin, Dimer
deg = 0.13; %Modify with degradation tag (for toxin only)
L = 1000; %Change as needed [uM]
options = odeset('RelTol', 1e-10);
tspan = [0 300];
y0 = [0; 0; 0];

% FIGURES HAD FONT SET TO SITKA TEXT AFTER CREATION
figure;
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("Concentrations of Holin, Antiholin, and Dimer: L = " + L/1000 + " mM");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
%ylim([0 100]);
axis auto
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");

% Bad code:
L = 10000;
figure;
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("Concentrations of Holin, Antiholin, and Dimer; L = " + L/1000 + " mM");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
%ylim([0 100]);
axis auto
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");

L = 30000;
figure;
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("Concentrations of Holin, Antiholin, and Dimer; L = " + L/1000 + " mM");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
%ylim([0 100]);
axis auto
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");

%%
% Analyze 3D
L_vec= linspace(0, 30, 15);
D_vec= logspace(-2, 0, 3);
output= zeros(length(L_vec), length(D_vec));
for i = 1:length(L_vec)
    for j= 1:length(D_vec)
        [t, x] = ode45(@(t,x) holinODE(t,x,L_vec(i) * 1000, D_vec(j)), tspan, y0, options);
        output(i,j) = getMolecule(x(end, 2));
        disp(i + ", " + j);
    end
end
disp(x)
figure;
imagesc(L_vec, D_vec, output);
xlabel("Lactate Concentration (mM)");
ylabel("Degradation constant (min^-^1)");
title("Number of Toxin Molecules");
colorbar;
caxis([0, 250])
%% with vs without degradation tag
deg = 0.13; %Modify with degradation tag (for toxin only)
L = 1000; %Change as needed [uM]
options = odeset('RelTol', 1e-10);
tspan = [0 300];
y0 = [0; 0; 0];

figure('DefaultAxesFontSize',14);
subplot(2,3,1)
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("L = " + L/1000 + " mM, d_{H} = 0.13 min^{-1}");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
ylim([0,1000])
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");
set(gca, 'FontName', 'PT Serif')

% Bad code:
L = 10000;
subplot(2,3,2)
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("L = " + L/1000 + " mM, d_{H} = 0.13 min^{-1}");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
ylim([0,1000])
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");
set(gca, 'FontName', 'PT Serif')

L = 30000;
subplot(2,3,3)
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("L = " + L/1000 + " mM, d_{H} = 0.13 min^{-1}");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
ylim([0,1000])
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");
set(gca, 'FontName', 'PT Serif')

deg = 0.0348;
subplot(2,3,4)
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("L = " + L/1000 + " mM, d_{H} = 0.0348 min^{-1}");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
ylim([0,1000])
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");
set(gca, 'FontName', 'PT Serif')

% Bad code:
L = 10000;
subplot(2,3,5)
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("L = " + L/1000 + " mM, d_{H} = 0.0348 min^{-1}");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
ylim([0,1000])
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");
set(gca, 'FontName', 'PT Serif')

L = 30000;
subplot(2,3,6)
[t, x] = ode45(@(t,x) holinODE(t,x,L,deg), tspan, y0);
disp(x);
plot(t, x, 'LineWidth', 1.5)
grid on
legend(["antiholin", "holin", "dimer"]);
legend('Location','best')
title("L = " + L/1000 + " mM, d_{H} = 0.0348 min^{-1}");
ylabel("Species concentration (μM)")
xlabel("Time (min)")
ylim([0,1000])
disp("[Lactate] = " + L/1000);
disp("Steady-state [Holin] = " + x(length(x),2) + " (μM)");
set(gca, 'FontName', 'PT Serif')

x0=10;
y0=10;
width=1200;
height=800;
set(gcf,'position',[x0,y0,width,height])
print(gcf,'with_without_deg_tag.png','-dpng','-r300');

%% 
function num_mol= getMolecule(uM)
    vol_E_coli= 1e-15;
    mole= 6.022e23;
    num_mol = (uM / 1e6) * (vol_E_coli * mole);
end
