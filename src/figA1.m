%% Figure A1: Solutions of the dispersion relation for K
% Luke Colosi | lcolosi@ucsd.edu | September 18th, 2021

%-------------------------------- Caption --------------------------------%
% 
%-------------------------------------------------------------------------%

clc, clear, close all;

% Set text interpreter 
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

% Set path for figures
fig_path = '../figs/'; 

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Set physical parameters
g = 9.81;                                                                   % Gravitational acceleration (units: m/s^2)
H = 4000;                                                                   % Water depth (units: m)
U = 0.55;                                                                   % Platform speed (units: m/s)
phi = 0;                                                                    % Platform propagation direction (units: degrees)

% Set wave parameters
f_ob = 0.5;                                                                 % Observed cyclicaly frequency (units: Hz)
omega = 2*pi*f_ob;                                                          % Observed radian frequency (units: rad/s)
k = -50:0.5:50;                                                             % x-component of wave number vector (units: rad/s)
l = -50:0.5:50;                                                             % y-component of wave number vector (units: rad/s)
[k_m, l_m] = meshgrid(k,l); 

% Compute wavenumber magnitude and direction
theta = atan2d(l_m,k_m);                                                    % Wave propagation direction (units: degrees)
K = sqrt(k_m.^2 + l_m.^2);                                                  % Wavenumber magnitude (units: rad/m)

% Compute the relative direction between the platform and waves
theta_r = theta - phi;                                                      % Angle between direction of platform and wave propagation (units: degrees)

% Compute dispersion relation in R reference frame (intrinsic frequency)
sigma_pos = (g*K.*tanh(K*H)).^(1/2); 
sigma_neg = -(g*K.*tanh(K*H)).^(1/2); 

% Compute dispersion relation in R' reference frame (observed frequency)
m = omega + K.*U.*cosd(theta_r); 

% Compute plane for threes relative heading cases
% plane_w = 
% plane_a = 
% plane_n = 

%% Plot dispersion relation surfaces 

clc, close all

% Create Figure and axes
fig = figure( 'Name', 'Dispersion relation surfaces');

% Set figure attributes
POS = [100 100 500 2000];                                                  
set(gcf,'color',[1 1 1])
set(gcf,'position',POS) 
fontsize = 12;
blue = [0 0.4470 0.7410]; 
red = [0.6350 0.0780 0.1840];
green = [0.4660 0.6740 0.1880];

% Plot surface
sc1 = surf(k,l,sigma_pos);

hold on 
    sc2 = surf(k,l,sigma_neg);
    sc3 = surf(k,l,m);
    % sc4 = surf(k,l,)
hold off

% Set figure attributes
colormap(flipud(cbrewer2('RdYlBu')))
shading flat
alpha(0.5)
sc1.FaceColor = blue;
sc2.FaceColor = red;
sc3.FaceColor = green;
xlabel('k (rad $\cdot$ m$^{-1}$)', 'Position', [23.991946187616122,-32.11125744837864,-37.56729252760999])
ylabel('l (rad $\cdot$ m$^{-1}$)', 'Position', [-41.00431899502064,22.52488724201885,-36.405558933024224])
zlabel('$\sigma$ or $\omega + \textbf{k}\cdot\textbf{U}$  (rad $\cdot$ s$^{-1}$)')
zlim([-40 40])
% set(gca,'ZScale','log')

