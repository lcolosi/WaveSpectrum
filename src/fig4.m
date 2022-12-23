%% Figure 4: Bifurcation frequency as a function of relative angle and platform speed
% Luke Colosi | lcolosi@ucsd.edu | August 17th, 2022

%-------------------------------- Caption --------------------------------%
% (a) Bifurcation frequency as a function of platform speed and the
% relative angle between platform heading and wave direction. Relative 
% angles range from 0 to 90 degrees correspond to the platform moving 
% with the waves. (b) Bifurcation frequency as a function of platform speed
% for a platform moving strictly in the direction of wave propagation.  
%-------------------------------------------------------------------------%

clc, clear, close all;

% Set default interpreter to latex for all text in figures
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

% Set path for figures
fig_path = '../figs/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
g = 9.81;                                                                   % Gravitational acceleration
u = 0:0.025:4;                                                              % Speed of wave glider (units: ms^-1) 
theta_r = 0:1:90;                                                           % Relative angle between mean current and wave direction (units: degrees)

% Initialize speed and relative angle parameters on a grid
[U, Theta_r] = meshgrid(u, theta_r);    

% Compute bifurcation frequency (units: Hz)
freq_b = (g)./(4*pi*U.*cosd(Theta_r));  

%% Plot bifurcation frequency in polar space 
close all

% Set variables for plotting: 
Ncircles = 5;                                                               % Number of azimuthal circles
pos = linspace(0, 4, Ncircles);                                             % Position of azimuthal circles and their labels
Contours = [0.1:0.1:2,5,10];                                                % Logarithmic scale colorbar tickmarks 
cmap = colormap(flipud(cbrewer2('RdBu', numel(Contours))));               % Colormap (flip so that red represents high f_b and blue represents low f_b) % , numel(levels)
idx_theta = theta_r == 0;                                                   % Index where theta vanishes
LineWidth = 0.5;
LineColor = 'k';
LineStyle = '-';
fontsize = 15;

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 1 1], 'Name', 'Birfurcation frequency in polar space');

%------------------- Subplot 1 -------------------%
subplot(1,2,1)
[~,cb] = polarcontourf(u, theta_r, freq_b, 'circlesPos', pos, 'Rscale',...
                       'linear', 'Ncircles', Ncircles, 'Nspokes', 10,...
                       'typeRose' ,'meteo', 'Contours', Contours,...
                       'colBar', 2, 'fontsize', fontsize, 'colormap', cmap,...
                       'LineWidth', LineWidth, 'LineColor', LineColor,...
                       'LineStyle', LineStyle);

% Set colorbar attributes 
cb.Ticks = log10([0.1, 0.2, 0.5, 1, 2, 5, 10]);
cb.TickLabels = {'0.1'; '0.2'; '0.5'; '1'; '2'; '5'; '10'};
cb.TickDirection = 'out';
cb.TickLength = 0.025;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = fontsize;
cb.Position = [0.4347    0.6399    0.0111    0.1985];
annotation('textbox',[0.4267 0.8434 0.0484 0.0374],'String',{'$f_c$ (Hz)'},...
           'Interpreter','latex','FontSize',15,'FontName','Helvetica Neue',...
           'EdgeColor','none');

% Set figure attributes
title('(a)')
set(gcf,'color',[1,1,1])
ax = gca;
ax.FontSize = fontsize;
set(gca,'TickLabelInterpreter','latex')
annotation('line', [0.13056,0.1222], [0.83434,0.83434])                     % Horizontal lines along vertical radial axis in subplot (a)
annotation('line', [0.13056,0.1222], [0.67503,0.67503])
annotation('line', [0.13056,0.1222], [0.51439,0.51439])
annotation('line', [0.13056,0.1222], [0.35374,0.35374])
annotation('line', [0.13125,0.1229], [0.19443,0.19443])
annotation('line', [0.12986,0.12986],[0.19579 0.1807])                      % Vertical lines along horizontal radial axis in subplot (a)
annotation('line', [0.21250,0.2125], [0.19712,0.1821])                      
annotation('line', [0.29514,0.2951], [0.19440,0.1794])
annotation('line', [0.37847,0.3785], [0.19445,0.1794])
annotation('line', [0.46042,0.4604], [0.19579,0.1807])
annotation('textbox',[0.3371 0.7577 0.0397 0.0375],...                      % Azimuthal axis tick label
           'String','$\theta_r \;(^\circ)$',...
           'Interpreter','latex','FontSize',15,...
           'FontName','Helvetica Neue','EdgeColor',[1 1 1]);
annotation('textbox',[0.10348,0.5019,0.0592,0.0375],...                     % Radial axis tick label
           'String',{'U (ms$^{-1}$)'},'Interpreter','latex',...
           'FontSize',15,'FontName','Helvetica Neue','EdgeColor',[1 1 1],...
           'Rotation', 90);
annotation('textbox',[0.1086 0.81923 0.0164 0.0281],...                     % Radial axis tick marks
           'String','4','Interpreter','latex','FontSize',15, ...
           'FontName','Helvetica Neue','EdgeColor',[1 1 1]);
annotation('textbox',[0.1107 0.6586 0.0164 0.0281],...
           'String','3','Interpreter','latex','FontSize',15,...
           'FontName','Helvetica Neue','EdgeColor',[1 1 1]);
annotation('textbox',[0.1093 0.5020 0.0164 0.0281],...
           'String','2','Interpreter','latex','FontSize',15,...
           'FontName','Helvetica Neue','EdgeColor',[1 1 1]);
annotation('textbox',[0.1100 0.3414 0.0164 0.0281],...
           'String','1','Interpreter','latex','FontSize',15,...
           'FontName','Helvetica Neue','EdgeColor',[1 1 1]);
annotation('textbox',[0.1107 0.1821 0.0164 0.0281],...
           'String','0','Interpreter','latex','FontSize',15,...
           'FontName','Helvetica Neue','EdgeColor',[1 1 1]);

% Set axis positions of figure
pos_ax1 = [0.1294, 0.1784, 0.4611 - 0.1294, 0.8303 - 0.1928]; % Note: position vector: [left bottom width height]

% Set variables for plotting annotations
vline = 0.025; % Length of vertical lines
line_s = 0.02; % Distance from figure 
linewidth = 1.5;

% Set position vector corresponding to the x-axis speed and y-axis
xpos = linspace(pos_ax1(1),pos_ax1(1)+pos_ax1(3),10000);
xspeed = linspace(0,4,10000);
ypos = pos_ax1(2)-line_s;

% Get xpos for the wave glider and saildrone
xpos_wg = xpos(xspeed >= 0 & xspeed <= 1);
xpos_sd = xpos(xspeed >= 0 & xspeed <= 3);

% Create platform speed annotations
%----- Wave Glider speed scale -----%

%-- Horizontal Line --%
xa = [xpos_wg(1) xpos_wg(end)];
ya = [ypos ypos];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)

%-- Vertical Lines --%
xa = [xpos_wg(1) xpos_wg(1)];
ya = [ypos-(vline/2) ypos+(vline/2)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)
xa = [xpos_wg(end) xpos_wg(end)];
ya = [ypos-(vline/2) ypos+(vline/2)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)

%-- Text box --%
str = '\bf{Wave Glider}';
dim = [0.1443 0.1480 0.0564 0.0295];
annotation('textbox',dim ,'String',str,'FitBoxToText','on','EdgeColor', 'none', 'Interpreter', 'latex');

%----- Sail Drone speed scale -----%

%-- Horizontal Line --%
xa = [xpos_sd(1) xpos_sd(end)];
ya = [ypos-(1.4*line_s) ypos-(1.4*line_s)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)

%-- Vertical Lines --%
xa = [xpos_sd(1) xpos_sd(1)];
ya = [ypos-(vline/2)-(1.4*line_s) ypos+(vline/2)-(1.4*line_s)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)
xa = [xpos_sd(end) xpos_sd(end)];
ya = [ypos-(vline/2)-(1.4*line_s) ypos+(vline/2)-(1.4*line_s)];
annotation('line',xa,ya,'Color','k', 'LineWidth', 1.5)

%-- Text box --%
str = '\bf{Saildrone}';
dim = [0.2353 0.1235 0.0497 0.0295];
annotation('textbox',dim ,'String',str,'FitBoxToText','on','EdgeColor', 'none', 'Interpreter', 'latex');

%------------------- Subplot 2 -------------------%
subplot(1,2,2)
semilogx(freq_b(idx_theta,:), u, '-k', 'LineWidth', 1)

% Set figure attributes
title('(b) $\theta_r = 0^\circ$')
xlabel('$f_c$ (Hz)')
ylabel('U (ms$^{-1}$)')
xlim([0.1,2])
grid on
set(gcf,'color',[1,1,1])
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', fontsize)
set(gca, 'position', [0.5090, 0.1928, 0.1382, 0.6466])

% Get position of figure
pos_ax2 = get(gca, 'Position'); 

% Set variables for plotting annotations
vline = 0.015; % Length of vertical lines
line_s = 0.02; % Distance from figure 

% Set position vector corresponding to the x-axis speed and y-axis
xpos = pos_ax2(1)+pos_ax2(3)+line_s;
yspeed = linspace(0,4,10000);
ypos = linspace(pos_ax2(2),pos_ax2(2)+pos_ax2(4),10000); 

% Get xpos for the wave glider and saildrone
ypos_wg = ypos(yspeed >= 0 & yspeed <= 1);
ypos_sd = ypos(yspeed >= 0 & yspeed <= 3);

% Create platform speed annotations
%----- Wave Glider speed scale -----%

%-- Horizontal Line --%
xa = [xpos xpos];
ya = [ypos_wg(1) ypos_wg(end)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)

%-- Vertical Lines --%
xa = [xpos-(vline/2) xpos+(vline/2)]; 
ya = [ypos_wg(1) ypos_wg(1)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)
xa = [xpos-(vline/2) xpos+(vline/2)]; 
ya = [ypos_wg(end) ypos_wg(end)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)

%-- Text box --%
str = '\bf{Wave Glider}';
dim = [0.6699 0.2284 0.0564 0.0295];
annotation('textbox',dim ,'String',str,'FitBoxToText','on','EdgeColor', 'none', 'Interpreter', 'latex', 'Rotation',90);

%----- Sail Drone speed scale -----%

%-- Horizontal Line --%
xa = [xpos+(1*line_s) xpos+(1*line_s)]; 
ya = [ypos_sd(1) ypos_sd(end)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)

%-- Vertical Lines --%
xa = [xpos-(vline/2)+(1*line_s) xpos+(vline/2)+(1*line_s)]; 
ya = [ypos_sd(1) ypos_sd(1)];
annotation('line',xa,ya,'Color','k', 'LineWidth', linewidth)
xa = [xpos-(vline/2)+(1*line_s) xpos+(vline/2)+(1*line_s)]; 
ya = [ypos_sd(end) ypos_sd(end)];
annotation('line',xa,ya,'Color','k', 'LineWidth', 1.5)

%-- Text box --%
str = '\bf{Saildrone}';
dim = [0.6895 0.3886 0.0453 0.0295];
annotation('textbox',dim ,'String',str,'FitBoxToText','on','EdgeColor', 'none', 'Interpreter', 'latex','Rotation',90);

% Save Figure
saveas(gcf, [fig_path 'figure_4.png'])