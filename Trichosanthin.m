%% Solving the ODEs
tspan = [0 20000]; %sec
y0 = 0;
[t, x] = ode45(@(t,x) therapeutic(t,x), tspan, y0);
t = t + 180; % time shift (export time)
x = x/3750; %adjust to cancer cell concentration

%% Plotting Figure
figure;
plot(t / (60 * 60), x, 'LineWidth', 1.5)
grid on
xlabel("Time (hours)") 
ylabel("Concentration (μM)")
title("Concentration of Trichosanthin in Cancer Cell")
legend("Trichosanthin")
legend("Location","best")

%% Steady State Concentrations
ss=x(length(t)); 
disp(ss)

%% Function script
function dRdt = therapeutic(t, x)

%parameters for dR/dt
basR = 0.387; %promoter basal expression(μM/s)
dR = 0.000385; %degradation rate of trich (s^-1)

dRdt = basR - dR * x;
end