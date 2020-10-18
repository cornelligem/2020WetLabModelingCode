%% define environment for model
 
%define viable radius of interest, includes tumor and some extracellular
%space
 
% since there are an average of 30 vessels/ mm^2, then 1/sqrt(3) mm on 
% average to a vessel
r = 25 + 1/sqrt(30); %CHANGE HERE https://cancerres.aacrjournals.org/content/63/17/5188.long
dispoints= 1000;
xmesh= linspace(0, r, dispoints);
 
tend= 1728000; % 1 week in seconds CHANGE HERE
tpoints= 1000;
tspan= linspace(0, tend, tpoints);
 
sol= pdepe(2, @pdefun, @icfun, @bcfun, xmesh , tspan);
 
%% plotting
 
%FIGURE FONT SET TO SITKA TEXT AFTER CREATION
figure(1); 
imagesc(xmesh / 10, tspan/ (60 * 60 * 24), sol * 1000); % CHANGE HERE this plot is not very readable
caxis([10^6 10^11])
hold on
a = colorbar;
a.Label.String = 'CFU/gram';
set(gca, 'ColorScale', 'log');
injection= 5 / 10; %your point goes here 
tumor= 25/10;
xline(injection, '--k', 'Linewidth', 1);
xline(tumor, '-k', 'Linewidth', 1);
xlabel("Position in tumor (cm)");
ylabel("Time after Injection (days)");
title('Bacterial Diffusion for Large Tumor (radius=25 mm)')
hold off
 
figure(2); 
imagesc(xmesh / 10, tspan/ (60 * 60 * 24), sol * 1000); % CHANGE HERE this plot is not very readable
hold on
a = colorbar;
a.Label.String = 'CFU/gram';
set(gca, 'ColorScale', 'log');
injection= 5 / 10; %your point goes here 
tumor= 25/10;
xline(injection, '--k', 'Linewidth', 1);
xline(tumor, '-k', 'Linewidth', 1);
xlabel("position in tumor (cm)");
ylabel("time after injection (days)");
xlim([tumor, r/10]);
hold off
 
% Animation of solution profile over time.
u_max = max(sol(:));  % determine maximal value of function u
figure(3);
plot(xmesh,sol(1,:),'b', xmesh, sol(end,:),'r');
axis([0 r 0 u_max])
title("Initial and End profiles");
xlabel('Coordinate x');
ylabel('u(x,t)');
 
%% Initial condition function
function u0= icfun(x)
    %sets initial concentration at 0
    %assume injection is in a spot of radius r and contains concentration c
    r= 5; %mm CHANGE THIS
    % mice injected with more than 2 million CFU suffered other side
    % effects, so keep it to 2 million CFU
    c= 2e6 / (r ^3); % CFU / mm^3
    if ( x > r)
        u0= 0;
    else 
        u0 = c;
    end
end
 
%% Boundary condition function
function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
    %assume injection is in a spot of radius r and contains concentration c
    r= 5; %CHANGE THIS
    c= 2e6 / (r ^3); % CFU / um^3
    % Set concentration of bacteria at boundaries to zero 
    pL = uL - c;
    qL = 0;
    pR = uR;
    qR = 0;
end
 
function [c, f, s] = pdefun(x, t, u, dudx)
    % define Mu based on x location, Mu is greater within the tumor than
    % outside of it
    r= 25; %radius of the tumor in mm CHANGE THIS
    lambda= 1; %look at growth curve CHANGE THIS
    A= 10^7; % CFU / mm^3 converted from 10^10 CFU/g https://www.nature.com/articles/s12276-019-0297-0
    % We will not account for change in tumor size because although the 
    % bacterial treatment does improve survival, papers suggest that tumors 
    % do not shrink significantly very quickly
    if x < r  %if inside the tumor, greater growth factor
        Mu= 1.736e-4; %CHANGE THIS
    else %else, slow growth
        Mu= 0; % CHANGE THIS
    end
    t_eff= invBacGrowthFun(Mu, lambda, A, u);
    dydt= bacGrowthFunDiff(Mu, lambda, A, t_eff);
    %define D
    D= 0.579e-6; % mm^2 CHANGE HERE
 
    %define differential equation as follows
    % (d/dt) u = D * (d^2 / dx^2) u + Mu * u
    c=1;
    f= D * dudx;
    s= dydt * u;
end
 
%% Get derivative of bacterial population given an effective time along the
% growth curve to determine how fast the bacteria will grow at a given
% concentration
% Equation found here: https://aem.asm.org/content/aem/56/6/1875.full.pdf
function dydt= bacGrowthFunDiff(mu, lambda, A, t)
    dydt= 4 * mu * exp( (4 * mu / A) * (lambda - t) + 2) / (1 + exp( (4 * mu / A) * (lambda - t) + 2))^2;
    if isnan(dydt)
        dydt= 0;
    end
end
 
%% get effective time in growth curve based on current concentration of
% bacteria, reverse logistic function with maximum density A, lag time
% lambda, and slope Mu. All the if statements are to deal with weird cases,
% like negative or NaN times. 
% Equation found here: https://aem.asm.org/content/aem/56/6/1875.full.pdf
function t= invBacGrowthFun(mu, lambda, A, y)
    if y==0
        t= 0;
    elseif (A / y - 1) <= 0
        t = Inf;
    else 
        t= lambda - (A / ( 4 * mu)) * (log(A / y - 1) - 2);
        if t < 0
            t = 0;
        end
    end
end
