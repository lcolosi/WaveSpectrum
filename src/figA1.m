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

%--------------------------------- Notes ---------------------------------%
% (1) Questions I need to answer and understand: 
%
%   (a) What is the mathematical justification for the k = 0 K solution? 
%       Answer: The + and minus signs are the result of a directional
%       ambiguity in the vanishing platform velocity case. 
%   (b) Why does the intrinsic frequency near the origin disappear in the 
%       l = 0 case or when phi is not equal to 0 or 180? 
%   (c) Why does the K solution for theta_r = 180 case have the right
%       magnitude but the wrong sign? 
%       Answer: The K solution has a sign ambiguity similar to k_p.
%       However, this ambiguity can be reloved by multiplying the
%       wavenumber magnitude by the cosine of relative theta. This provides
%       us with the wavenumber components only the l = 0 plane as shown in
%       panels b and c. 
%
% (2) Improvements to code: 
% 
%   (a) Consider other ways to illustrate the K solutions with varying U or
%       theta_r in the intrinsic or observed frequency space. Idea: create
%       a figure of f_ob as a function of K with various theta_r. 
%           DONE! (see development code)
% 
%   (b) Create a figure of group and phase speed as a function of intrinsic
%       frequency with a horizontal line for the uniform platform speed. 
%           DONE! (see development code)
% 
%   (c) Automate code so the proper projection projection onto the k = 0
%       and l =0 plane for the k and l components are obtained. 
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
legend([sc1 sc2, sc3], {'$\sigma(\textbf{k})$'; '-$\sigma(\textbf{k})$';'$\omega + \textbf{k}\cdot\textbf{U}$'}, 'Position', [0.390457885867418,0.448655914170768,0.054062316152785,0.050268817012028], 'FontSize', fontsize-4) 
set(gca,'FontSize',fontsize)

% Set for direction of platform propagation
annotation('arrow',[0.296527777777778 0.338888888888889],[0.188516129032258 0.204301075268817]);
annotation('textbox',[0.304472222222222 0.201612903225807 0.018875365787082 0.029569892473118],...
    'String',{'$\textbf{U}$'},'Interpreter','latex', ...
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
legend([sc1 sc2, sc3], {'$\sigma(\textbf{k})$'; '-$\sigma(\textbf{k})$';'$\omega + \textbf{k}\cdot\textbf{U}$'}, 'Position', [0.662377002938179 0.605141248535899 0.0626229970618208 0.0633333331795148], 'FontSize', fontsize-4) 
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
saveas(gcf, [fig_path 'figA1.png'])

%% Development code

%------------------- Observed Frequency Figure -------------------%

% Compute the observed frequency 
omega_pos = sigma_pos + K.*U.*cosd(theta_r);

%------------ Set plotting variables ------------%
% Obtain RGB triplet for colormap
cmap = colormap(cbrewer2('RdYlBu'));

% Set relative angle for plotting the number of relative angles plotted
theta_r_p = 0:30:180;

% Bin average colorbar to match number of relative angle cases
ntriplet = floor(length(cmap)/length(theta_r_p))*length(theta_r_p);             % number of RGB triplets that divide evenly into number of angle bins
ncolor = floor(length(cmap)/length(theta_r_p));                               % number of triplet to be bin averaged
cmap_bin = sepblockfun(cmap(1:ntriplet,:), [ncolor,1], @mean);

% Set colormap and colorbar variables
z_m = linspace(0,180,7);
z = 0:12.8571:180;

% Obtain RGB triplet of colormap
cmap = colormap(cmap_bin);

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 0.5 1]);

% Loop through discrete angles
for iangle = 0:30:180

    % Set the indicies for the wave encounter angles
    idx = round(theta_r,0) == iangle;

    % Set the indicies colorbar
    idx_color = theta_r_p == iangle;

    % Set line color
    [~,idx_c] = min(abs(z_m-theta_r_p(idx_color)));

    % Plot observed frequency trace at ith wave encounter angle
    hold on 
    plot(K(idx),omega_pos(idx), '-', 'LineWidth', 2, 'Color', cmap(idx_c,:));
    hold off

    % Plot intrinsic frequency curve
    if iangle == 180
        hold on
        plot(K(idx),sigma_pos(idx), '--', 'LineWidth', 2, 'Color', 'k');
        hold off
    end

end

% Set figure attributes
axis equal 
axis tight
xlabel('$K$ (rad m$^{-1}$)','Interpreter','Latex')
ylabel('$\omega$ (rad s$^{-1}$)','Interpreter','Latex')
xlim([0 50]); ylim([-10 50])
xticks(0:10:50); yticks(-10:10:50);
grid on
set(gca,'FontSize',fontsize)
set(gca,'Box','on')

%-------- Set Colorbar --------%
colormap(cmap_bin);
cb = colorbar; 
cb.Ticks = z(2:2:end)/max(z);
cb.TickLabelInterpreter = 'Latex'; 
cb.Label.Interpreter = 'Latex';
cb.TickLabels = {'0', '30', '60', '90', '120', '150', '180'};   
cb.FontSize = fontsize;
%cb.TickLabels = num2str(theta_r');
cb.Label.String = '';
cb.TickDirection = 'out';
annotation('textbox',[0.838694444444445 0.9345311827957 0.108 0.0359000000000002],...
           'String',{'$\theta_r \;(^\circ)$'},'Interpreter','latex',...
           'FontSize',fontsize,'EdgeColor','none');


%------------------- Group and phase speed Figure -------------------%

% Set index for intrinsic frequency in the theta_r = 0 direction
idx_zero = theta_r == 0; 

% Compute ground and phase velocity
c_g = (1/2)*(g./sigma_pos(idx_zero)); 
c_p = (g./sigma_pos(idx_zero)); 

% Create Figure and axes
fig = figure('units','normalized','outerposition',[0 0 0.5 1]);

% Plot group and phase velocity
pc1 = plot(sigma_pos(idx_zero), c_g, '-', 'LineWidth', 2, 'color', blue);
hold on 
pc2 = plot(sigma_pos(idx_zero), c_p, '-', 'LineWidth', 2, 'color', red);
pc3 = yline(U, '-', ['$U =$ ' num2str(U) ' ms$^{-1}$'], 'LineWidth', 2, 'color', 'k', 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'FontSize',fontsize);
hold off

% Set figure attributes

xlabel('$\sigma$ (rad s$^{-1}$)','Interpreter','Latex')
ylabel('(m s$^{-1}$)','Interpreter','Latex')
grid on
set(gca,'FontSize',fontsize)
set(gca,'Box','on')
leg = legend([pc1, pc2], '$c_g$', '$c_p$','location', 'northwest', 'fontsize', fontsize-3);

