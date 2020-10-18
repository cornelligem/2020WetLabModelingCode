%% Solving the ODEs
tspan =[0 200]; 
c0=0; %initially 0 μM mCardinal protein
[t,y] = ode45(@(t,y) mCard(t,y), tspan, c0); 

%% Plotting Figure
figure
plot(t,y, 'LineWidth',1.5) 
grid on
xlabel("Time (min)")
ylabel("Concentration (μM)")
title("mCardinal Concentration over Time")
legend("mCardinal")
legend("Location","best")

%% Steady State Concentrations
ss=y(length(t)); 
disp(ss)

%% function for mCardinal modeling
function dxdt= mCard(t, x)

%parameters for dC/dt
basE = 23.22; %basal expression (μM/min)
delC= 0.0575; %mCardinal degradation rate (min^-1)

dxdt = basE - delC*x;

end