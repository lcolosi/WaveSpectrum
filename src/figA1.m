%% Figure A1: Solutions of the dispersion relation for K
% Luke Colosi | lcolosi@ucsd.edu | September 18th, 2021

%-------------------------------- Caption --------------------------------%
% (a) An illustration of the surface of revolution \ref{dispersion_surface}
% (blue ($+\sigma$) and red ($-\sigma$) surfaces) and the plane
% \ref{doppler_velocity_plane} (green surface) with $\omega = 3.1412$
% rad s$^{-1}$, $U = 0.55$ ms$^{-1}$, and $\phi = 0^\circ$. (b) and (c)
% Values of \ref{dispersion_surface} and \ref{doppler_velocity_plane} on
% the surfaces of constant wavenumber for $l= 0$ and $k = 0$, respectively.
% $k$ solutions of the dispersion relation, using variables defined in
% Table~\ref{t1}, are at intersection points of \ref{dispersion_surface}
% and \ref{doppler_velocity_plane} (shown as the vertical black lines).
% The platform's velocity vector $U$ is shown in the \textbf{k}-plane. 
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
phi = 0;                                                                    % Platform propagation direction (units: degrees; Directional convention: CCW, reference east, going towards)

% Set wave parameters
f_ob = 0.5;                                                                 % Observed cyclical frequency (units: Hz)
omega = 2*pi*f_ob;                                                          % Observed radian frequency (units: rad/s)
k = -50:0.5:50;                                                             % x-component of wave number vector (units: rad/s)
l = -50:0.5:50;                                                             % y-component of wave number vector (units: rad/s)
[k_m, l_m] = meshgrid(k,l); 

% Compute wavenumber magnitude and direction
theta = atan2d(l_m,k_m);                                                    % Wave propagation direction (units: degrees; Directional convention: CCW, reference east, going towards)
K = sqrt(k_m.^2 + l_m.^2);                                                  % Wavenumber magnitude (units: rad/m)

% Convert directional convention to CW, reference north, going towards
theta = mod(90 - theta,360);                                     
phi = mod(90 - phi, 360);

% Compute the relative direction between the platform and waves
theta_r = mod(theta - phi,360);                                             % Angle between direction of platform and wave propagation (units: degrees)

% Compute dispersion relation in R reference frame (intrinsic frequency)
sigma_pos = (g*K.*tanh(K*H)).^(1/2); 
sigma_neg = -(g*K.*tanh(K*H)).^(1/2); 

% Compute dispersion relation in R' reference frame (observed frequency)
m = omega + K.*U.*cosd(theta_r);

% Compute the l = 0 surface
%l_zero = 

% Compute analytic K solutions for intersection point along k = 0 and l = 0
% planes (pick the k = 0 and l = 0 relative angles from theta_r). 
theta_r_kp = [0, 180]; 
theta_r_lp = [90,270]; 

%-------------------- l = 0 plane --------------------%
K_lt = (g - 4*pi*f_ob*U*cosd(theta_r_kp(1)) - sqrt(g^2 - 8*pi*f_ob*g*U*cosd(theta_r_kp(1))))/(2 * U^2 * cosd(theta_r_kp(1))^2);
K_st = (g - 4*pi*f_ob*U*cosd(theta_r_kp(1)) + sqrt(g^2 - 8*pi*f_ob*g*U*cosd(theta_r_kp(1))))/(2 * U^2 * cosd(theta_r_kp(1))^2);
K_r = (g - 4*pi*(-f_ob)*U*cosd(theta_r_kp(1)) + sqrt(g^2 - 8*pi*(-f_ob)*g*U*cosd(theta_r_kp(1))))/(2 * U^2 * cosd(theta_r_kp(1))^2);
K_a = (g - 4*pi*f_ob*U*cosd(theta_r_kp(2)) - sqrt(g^2 - 8*pi*f_ob*g*U*cosd(theta_r_kp(2))))/(2 * U^2 * cosd(theta_r_kp(2))^2);

%-------------------- k = 0 plane --------------------%
Kp = ((2*pi*f_ob)^2)/g;

% Compute the wavenumber components k and l in the k = 0 and l = 0 planes

%-------------------- l = 0 plane --------------------%
k_lt = cosd(theta_r_kp(1))* K_lt;
k_st = cosd(theta_r_kp(1))* K_st;
k_r = cosd(theta_r_kp(2))* K_r;
k_a = cosd(theta_r_kp(2))* K_a;

%-------------------- k = 0 plane --------------------%
lp_p = sind(theta_r_lp(1))*Kp; 
ln_p = sind(theta_r_lp(2))*Kp; 

% Set the indicies for the wave encounter angles
ind_parallel = theta_r == 0 | theta_r == 180;  
ind_normal = theta_r == 270 | theta_r == 90;

%% Plot dispersion relation surfaces 
clc, close all; 

% Create Figure and axes
fig = figure('units','normalized','outerposition',[0 0 1 1]);

% Set plotting variables
blue = [0 0.4470 0.7410]; 
red = [0.6350 0.0780 0.1840];
green = [0.4660 0.6740 0.1880];
black = [0,0,0];
fontsize = 21;
set(gcf,'color',[1 1 1])

%------------------- Dispersion Relation Surfaces -------------------%
ax1 = subplot(2,2,[1,3]);

% Plot dispersion relation surfaces
sc1 = surf(k,l,sigma_pos);

hold on 
    sc2 = surf(k,l,sigma_neg);
    sc3 = surf(k,l,m);
hold off

% Set figure attributes
title('(a)', 'Position', [13.907391675107192,2.577847801991084,49.074404163666486])
shading flat
alpha(0.5)
sc1.FaceColor = blue; sc2.FaceColor = red; sc3.FaceColor = green;
xlabel('$k$ (rad m$^{-1}$)', 'Position', [13.57557284106133,-27.52171618968987,-50.8773940577874])
ylabel('$l$ (rad m$^{-1}$)', 'Position', [-34.02161915211582,14.058365548812617,-50.99054550473954])
zlabel('m (rad s$^{-1}$)')
xlim([-50 50]); ylim([-50 50]); zlim([-40 40])
xticks(-50:25:50); yticks(-50:25:50); zticks(-40:10:40)
view([-36.8410110594965 9.62376237623763]);
legend([sc1 sc2, sc3], {'$\sigma(\textbf{k})$'; '-$\sigma(\textbf{k})$';'$\omega + \textbf{k}\cdot\textbf{u}$'}, 'Position', [0.390457885867418,0.448655914170768,0.054062316152785,0.050268817012028], 'FontSize', fontsize-4) 
set(gca,'FontSize',fontsize)

% Set for direction of platform propagation
annotation('arrow',[0.296527777777778 0.338888888888889],[0.188516129032258 0.204301075268817]);
annotation('textbox',[0.304472222222222 0.201612903225807 0.018875365787082 0.029569892473118],...
    'String',{'$\textbf{u}$'},'Interpreter','latex', ...
    'HorizontalAlignment','center','EdgeColor','none', 'FontSize', fontsize);

%------------------- Dispersion Relation curves: theta_r = 0, 180 or l = 0 plane -------------------%
ax2 = subplot(2,2,2);

% Plot dispersion relation curve
sc1 = plot(k_m(ind_parallel),sigma_pos(ind_parallel), 'color', blue, 'LineWidth',1.5);

hold on 
    sc2 = plot(k_m(ind_parallel),sigma_neg(ind_parallel), 'color', red, 'LineWidth',1.5);
    sc3 = plot(k_m(ind_parallel),m(ind_parallel), 'color', green, 'LineWidth',1.5);
    xline(k_lt, '-k', '$k_{lt}$', 'LineWidth', 1, 'Interpreter','latex','LabelOrientation','aligned', 'LabelHorizontalAlignment','right', 'LabelVerticalAlignment', 'top', 'FontSize',fontsize, 'Color', black)
    xline(k_st, '-k', '$k_{st}$', 'LineWidth', 1, 'Interpreter','latex','LabelOrientation','aligned', 'LabelHorizontalAlignment','right', 'LabelVerticalAlignment', 'top', 'FontSize',fontsize, 'Color', black)    
    xline(k_a, '-k', '$k_a$', 'LineWidth', 1, 'Interpreter','latex','LabelOrientation','aligned', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment', 'top', 'FontSize',fontsize, 'Color', black)
    xline(k_r, '-k', '$k_r$', 'LineWidth', 1, 'Interpreter','latex','LabelOrientation','aligned', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment', 'top', 'FontSize',fontsize, 'Color', black)
hold off

% Set figure attributes
title('(b) $\theta_r = 0^\circ, 180^\circ$')
xlabel('$k$ (rad m$^{-1}$)')
ylabel('m (rad s$^{-1}$)')
xlim([-50 50]); ylim([-40 40])
xticks(-50:25:50); yticks(-40:20:40);
legend([sc1 sc2, sc3], {'$\sigma(\textbf{k})$'; '-$\sigma(\textbf{k})$';'$\omega + \textbf{k}\cdot\textbf{u}$'}, 'Position', [0.662377002938179 0.605141248535899 0.0626229970618208 0.0633333331795148], 'FontSize', fontsize-4) 
grid on
set(gca,'FontSize',fontsize)

%------------------- Dispersion Relation curves: theta_r = 90, 270 or k = 0 plane -------------------%
ax3 = subplot(2,2,4);

% Plot dispersion relation curve
sc1 = plot(k_m(ind_parallel),sigma_pos(ind_parallel), 'color', blue, 'LineWidth',1.5);

hold on 
    sc2 = plot(k_m(ind_parallel),sigma_neg(ind_parallel), 'color', red, 'LineWidth',1.5);
    sc3 = plot(l_m(ind_normal),m(ind_normal), 'color', green, 'LineWidth',1.5);
    xline(ln_p, '-k', '$k_p^{-}$', 'LineWidth', 1, 'Interpreter','latex','LabelOrientation','aligned', 'LabelHorizontalAlignment','left', 'FontSize',fontsize, 'Color', black)
    xline(lp_p, '-k', '$k_p^{+}$', 'LineWidth', 1, 'Interpreter','latex','LabelOrientation','aligned', 'LabelHorizontalAlignment','right', 'FontSize',fontsize, 'Color', black)    
    
hold off

% Set figure attributes
title('(c) $\theta_r = 90^\circ, 270^\circ$')
xlabel('$l$ (rad m$^{-1}$)')
ylabel('m (rad s$^{-1}$)')
xlim([-50 50]); ylim([-40 40])
xticks(-50:25:50); yticks(-40:20:40);
grid on
set(gca,'FontSize',fontsize)

% Save Figure
print(fig,'-depsc', [fig_path 'figA1.eps'], '-r300');                        % Resolution units: dpi (dots per inch); default resolution: 150 dpi 
