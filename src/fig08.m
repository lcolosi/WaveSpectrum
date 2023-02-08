%% Figure 9:  Observations of Wave Spectra power spectral density modulation 
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% (a)  Wave Glider Stokes' heading (using the coming from directional
% convention) for the small 500 m box legs. Red triangles denote the legs
% where the platform is moving in the direction of wave propagation for
% high-frequency waves ($f_{ob} > 10^{-1}$ Hz). Blue triangles denote the
% legs where the platform is moving against the direction of wave
% propagation for high-frequency waves. These high-frequency waves have
% a mean wave direction coming from the Northwest (approximately
% $300^\circ$). (b) Saturation spectrogram for Wave Glider Stokes'
% repeated small 500 m box formation with red and blue colored ovals
% along the top axis specifying the relative angle between platform
% heading and mean wave direction for high-frequency wave using the
% same color code as in Figure~\ref{f9}a. Saturation spectra within the
% two time periods outlined by dot-dashed vertical lines correspond to
% the averaged spectra in Figure~\ref{f9}c. (c) Two averaged saturation
% spectra showing against (red curve) and with (blue curve) wave cases.
%-------------------------------------------------------------------------% 

clc, clear, close all;

% Set text interpreter 
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

% Set path to data
%%%% Vehicle %%%%.
vehicle = 'STOKES';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';

% Set path for figures
fig_path = '../figs/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
g = 9.81;                                                                   % gravitational acceleration (m/s) 
dir_con = {'CW', 'cf', 'rn'};                                               % Directional conventions

% Temporal parameters
period = 10;                                                                % Time duration for computing wave spectra
date_i = '09-Sep-2020 02:30:00';                                            % Start date 
date_f = '11-Sep-2020 16:10:00';                                            % End date
date_o = {date_i, date_f}; 
%------- Novatel ------- %
fe_n = 20;                                                                  % Sampling rate (Hz)
dt_n = 1/fe_n;                                                              % Time interval between measurements (s)
%------- Weather Station ------- %
fe_w = 1; 
dt_w = 1/fe_w; 

% Spectral parameters                                                      
fn = (1/2)*fe_n;                                                            % Nyquist Frequency (Hz)
df = 0.01;                                                                  % Frequency resolution
nfft = fe_n/df;                                                             % Window length 
f = (0:df:fn);                                                              % Observed frequency (Hz)
dtheta = 5;                                                                 % Angular resolution (degrees)
ntheta = 360/dtheta+1;                                                      % Number of angles
lambda_c = 1.5613;                                                          % Wavelength cutoff (meters)
f_noise = sqrt(g/(2*pi*lambda_c));                                          % Noise frequency cutoff 
toolbox = 'WAFO';                                                           % Method used to compute directional spectrum 
variables = 'heave_velocity';                                               % Heave and horizontal velocity are used to compute the direction spectrum
scaling = false;                                                            % Variance of directional spectrum is not scaled to match variance of heave spectrum 

% Upload and process novatel and weather station data
[nov_s, ~, nlegs_s] = process_wg_data(vehicle, ROOT, date_o, dir_con,...
                                     [dt_n, dt_w]);

%% Compute Wave Spectra

% Loop through legs
for n = 1:nlegs_s
    
    % Compute Directional Spectrum  
    eval(['[nov_s.Sd(:,:,n), nov_s.f, nov_s.theta] = compute_directional_spectrum(nov_s.L' num2str(n) '.heave, nov_s.L' num2str(n) '.VEL_east, nov_s.L' num2str(n) '.VEL_north, nov_s.L' num2str(n) '.VEL_up, nov_s.L' num2str(n) '.time_20hz, f, df, dtheta, nfft, fe_n, toolbox, variables, scaling, dir_con);']) 
    
    % Compute the Omni-directional Spectra from directional wave spectra 
    nov_s.spectrogram_omni(:,n) = sum(nov_s.Sd(:,:,n) * dtheta, 1);

    % Compute Mean and standard deviation for direction of platform
    eval(['[nov_s.mD_legs(n), nov_s.stdD_legs(n)] = direction_stats(nov_s.L' num2str(n) '.true_course, dt_n, 0);']) % Convention: CW, coming from, ref north
 
    % Compute the half way point of the time series
    eval(['nov_s.time_legs(n) = nov_s.L' num2str(n) '.time_20hz(round(length(nov_s.L' num2str(n) '.time_20hz) - (1/2)*length(nov_s.L' num2str(n) '.time_20hz)));']) 
    
    % Compute mean latitude
    eval(['nov_s.mlat_legs(n) = mean(nov_s.L' num2str(n) '.latitude);'])

end

% Obtain frequency indicies below high frequency cut off
Inoise_s = find(nov_s.f <= f_noise);

% Remove frequencies and beyond 1 Hz 
nov_s.f_ob = nov_s.f(Inoise_s);
nov_s.spectrogram_omni_f_ob = nov_s.spectrogram_omni(Inoise_s,:);
nov_s.Sd_f_ob = nov_s.Sd(:,Inoise_s,:);

% Compute equivalent saturation spectrum 
nov_s.sat_spectrogram_omni_f_ob = nov_s.spectrogram_omni_f_ob .* (nov_s.f_ob').^(5);

%% Grab data from small box trajectory and split into with and against waves legs

% Obtain latitude indicies for small box cutoff
Ilat_s = find(nov_s.mlat_legs < 32.93);  % Indices for legs in small box

% Grab mean course direction, time and saturation spectrogram for small box
mD_small = nov_s.mD_legs(Ilat_s); 
time_small = nov_s.time_legs(Ilat_s);
sat_spec_small = nov_s.sat_spectrogram_omni_f_ob(:,Ilat_s);

% Make small adjustment to mean direction for plotting 
nov_s.mD_legs(nov_s.mD_legs > 350) = 0;                                     % Set mean directions near 360 to 0 (directions wrap)
mD_small(mD_small > 350) = 0;                                               

% Split mean direction time and spectrogram into directions with and into
% waves (!!!USING COMING FROM CONVENTION!!!)

% Initialize wave parameters 
mwd = 300;                                                                  % Mean Wave Direction for high frequency wave system (directional convention set by dir_con)
Edges = mod([270,90] - (360 - mwd),360);                                    % Boundaries for with and against waves azimuthal regimes 

%------ With Waves ------%
idx_with = mD_small > Edges(1) | mD_small < Edges(2);
mD_with = mD_small(idx_with); 
time_with = time_small(idx_with);
sat_spec_with = sat_spec_small(:,idx_with);
%------ Against Waves ------%
idx_against = mD_small > Edges(2) & mD_small < Edges(1);
mD_against = mD_small(idx_against); 
time_against = time_small(idx_against);
sat_spec_against = sat_spec_small(:,idx_against);

% Compute average saturation spectra for a given box
av_sat_spec_with = mean(sat_spec_small(:,8:9),2);
av_sat_spec_against = mean(sat_spec_small(:,6:7),2);

%% Plot mean platform direction saturation wave spectrogram, and case study of two saturation spectrum 

% Set plotting variables
a = 1:2:length(Ilat_s); 
b = 0:2:length(Ilat_s);
idx_with = sort([1,a(3:2:length(a)),b(3:2:length(b))]); 
idx_into = sort([a(2:2:length(a)),b(2:2:length(b))]);
first = datetime(2020, 09, 10,23,0,0);
last = datetime(2020, 09, 11,13,0,0);
t_ticks = datenum(first:hours(1):last);
t_case = find(nov_s.time_legs > datenum(2020, 09, 11,2,0,0) & nov_s.time_legs < datenum(2020, 09, 11,3,20,0));
fontsize = 19;
red = [0.6350 0.0780 0.1840]; 
blue = [0 0.4470 0.7410];

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 1 0.7])

%------------- Subplot 1 -------------%
ax1 = subplot(2,2,1);

% Plot mean direction per leg
plot(nov_s.time_legs(Ilat_s), nov_s.mD_legs(Ilat_s), '-k', 'LineWidth', 2);

hold on 
    pc1 = plot(time_against, mD_against, '^', 'MarkerSize', 10, 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k'); 
    pc2 = plot(time_with, mD_with, '^', 'MarkerSize', 10, 'MarkerFaceColor', red, 'MarkerEdgeColor', 'k');
hold off

% Set figure attributes
title('(a)')
ylabel('Wave Glider Heading ($^{\circ}$)', 'Interpreter', 'latex')
xlim([nov_s.time_legs(Ilat_s(1)), nov_s.time_legs(end)])
ylim([0, 360])
yticks(linspace(0,360,5))
xticks(datenum(t_ticks(2:end)))
datetick('x', 'HH', 'keepticks')
grid on 
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
l = legend([pc1, pc2], {'Against Waves', 'With Waves'}, 'Interpreter', 'latex', 'Location', 'best', 'fontsize', 14);
rect = [0.334161689302536,0.849899970974651,0.091234673394097,0.069147004406162];
set(l, 'Position', rect)

%------------- Subplot 2 -------------%
ax2 = subplot(2,2,3);

% Plot Spectrogram
pc = pcolor(time_small, nov_s.f_ob, sat_spec_small);

hold on

% Loop through legs 
for ileg = 1:length(Ilat_s)-4
    
    % Plot rectangle above leg
    if sum(ileg == idx_with)
        rectangle('Position',[time_small(ileg) 1 time_small(ileg+1) - time_small(ileg) 0.2], 'FaceColor', red, 'EdgeColor', 'k', 'Clipping', 'off', 'Curvature', 1)
    elseif sum(ileg == idx_into) 
        rectangle('Position',[time_small(ileg) 1 time_small(ileg+1) - time_small(ileg) 0.2], 'FaceColor', blue, 'EdgeColor', 'k', 'Clipping', 'off', 'Curvature', 1)
    end
    
end

xline(time_small(6), '-.', 'LineWidth', 2, 'Color', 'k', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', 15)
xline(time_small(8), '-.', 'LineWidth', 2, 'Color', 'k', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', 15)
xline(time_small(10), '-.', 'LineWidth', 2, 'Color', 'k')

hold off

% Set figure attributes
title('(b)')
xlabel('UTC time from 11 Sep 2020 (hrs)')
ylabel('$f_{ob}$ (Hz)')
xlim([nov_s.time_legs(Ilat_s(1)), nov_s.time_legs(end)])
ylim([3*10^-2 10^0])
xticks(datenum(t_ticks(2:end)))
datetick('x', 'HH', 'keepticks')
set(get(gca,'title'),'Position',[nov_s.time_legs(61) 1.2 0.10011])
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
pc.EdgeColor = 'none';
%annotation('rectangle',[0.172916666666667 0.116152450090744 0.015277777777778 0.334509240835089],...
%           'LineStyle','none','FaceColor',blue,'FaceAlpha',0.3);
%annotation('rectangle',[0.189657407407406 0.114337568058076 0.011037037037039 0.336800316814746],...
%           'LineStyle','none','FaceColor',red,'FaceAlpha',0.3);

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = 'B$_{ob}$($t,f_{ob}$) (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

%------------- Subplot 3 -------------%
subplot(2,2,[2,4]);

% Plot average saturation spectra for one box trajectory  
pc1 = loglog(nov_s.f_ob, av_sat_spec_against, '-', 'LineWidth', 2, 'color', blue);

hold on
pc2 = loglog(nov_s.f_ob, av_sat_spec_with, '-', 'LineWidth', 2, 'color', red);
hold off

% Set figure Attributes
title('(c)')
ylabel('B$_{ob}$($f_{ob}$) (m$^2$ Hz$^{4}$)','Interpreter','Latex')
xlabel('$f_{ob}$ (Hz)','Interpreter','Latex')
xlim([3*10^-2 10^0])
grid on 
axis square
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'TickLabelInterpreter','latex')
set(gca,'Box','on')
set(gca, 'color', 'w')

%-------- Set the position of the subplots --------%
% find current position [x,y,width,height]
pos1 = get(ax1,'Position');
pos2 = get(ax2,'Position');

% set width of second axes equal to first
pos1(3) = pos2(3);
set(ax1,'Position',pos1)

% Save Figure
saveas(gcf, [fig_path 'fig09.png'])