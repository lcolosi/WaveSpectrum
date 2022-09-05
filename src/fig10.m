%% Figure 10: Full time series saturation spectrograms for Stokes DELMAR2020 and WHOI43 SMODE2021 Experiments
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% Observed omni-directional (a) and saturation (b) wave spectrograms 
% for the Stokes wave glider.  
%-------------------------------------------------------------------------% 

clc, clearvars -except DM SM, close all;

% Set text interpreter 
set(0,'defaultTextInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

% Set path to data
%%%% Vehicle %%%%.
vehicle1 = 'STOKES';
vehicle2 = 'WHOI43';

%%%% Root %%%%
ROOT_dm = '../data/DELMAR2020/';
ROOT_sm = '../data/SMODE2021/';

% Set path for figures
fig_path = '../figs/';

%% Call Data 

%--------- DELMAR 2020 ---------% 
if isempty(whos('DM')) 

    DM = load([ROOT_dm vehicle1 '/DELMAR2020_spec_' vehicle1 '.mat']); % Mapping variables computed in delmar_intrinsic_freq_spectrogram.m file

end

%--------- SMODE Pilot 2021 ---------% 
if isempty(whos('SM')) 

    SM = load([ROOT_sm vehicle2 '/SMODE2021_spec_' vehicle2 '.mat']); % Mapping variables computed in smode_intrinsic_freq_spectrogram.m file

end


%% Plot saturation wave spectrogram before and after mapping 
clc, close all; 

% Set power spectral densities equal to zero to NaN
%------- DM -------%
idx_zero_dm = DM.spec_ob == 0; 
DM.spec_ob(idx_zero_dm) = NaN; DM.spec_in(idx_zero_dm) = NaN; 

%------- SM -------%
idx_zero_sm = SM.spec_ob == 0; 
SM.spec_ob(idx_zero_sm) = NaN; SM.spec_in(idx_zero_sm) = NaN; 

% Set frequency low cutoff
Ifreq = find(DM.f > 0.03);                                                  %Frequency low cutoff 

% Set plotting parameters
t_ticks_dm = datetime('09-Sep-2020 12:00:00'):hours(12):datetime('11-Sep-2020 12:00:00');
t_ticks_sm = datetime('29-Oct-2021 00:00:00'):days(1):datetime('04-Nov-2021 00:00:00'); 
fontsize = 16;

% Create figure
figure('units','normalized','outerposition',[0 0 1 0.6])

%------------ Subplot 1 ------------% 
subplot(2,2,1)

% Plot observed saturation spectrogram from DELMAR2020
pc = pcolor(DM.time, DM.f(Ifreq), DM.spec_ob(Ifreq,:));

% Set figure attributes
title('(a)')
pc.EdgeColor = 'none';
ylabel('$f_{ob}$ (Hz)', 'Interpreter', 'latex')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca, 'TickLength', [0.007, 0.007]) 
xticks(datenum(t_ticks_dm))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([DM.time(2), DM.time(end)])
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca,'Box','on')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu', 100)))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{ob}^5 \cdot$ S($t,f_{ob}$) (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [ 10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ;  
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

%------------ Subplot 2 ------------% 
subplot(2,2,2)

% Plot observed saturation spectrogram from SMODE2021
pc = pcolor(SM.time, SM.f(Ifreq), SM.spec_ob(Ifreq,:));

% Set figure attributes
title('(b)')
pc.EdgeColor = 'none';
xlim([SM.time(1), SM.time(end)])
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca, 'TickLength', [0.007, 0.007]) 
xticks(datenum(t_ticks_sm))
datetick('x', 'mmm dd', 'keepticks')
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca,'Box','on')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu', 50)))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{ob}^5 \cdot$ S($t,f_{ob}$) (m$^2$ Hz$^{4}$)';
caxis([10^-7, 2*10^-3]);
cb.Ticks = [10^-7; 10^-6; 10^-5; 10^-4; 10^-3];
cb.TickLabels = {'$10^{-7}$'; '$10^{-6}$'; '$10^{-5}$'; '$10^{-4}$'; '$10^{-3}$'}; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

%------------ Subplot 3 ------------% 
subplot(2,2,3)

% Plot intrinsic frequency saturation spectrogram from DELMAR2020
pc = pcolor(DM.time, DM.f(Ifreq), DM.spec_in(Ifreq,:));

% Set figure attributes
title('(c)')
pc.EdgeColor = 'none';
xlabel('UTC time from Sep 9$^{\textrm{th}}$, 2020', 'Interpreter', 'latex')
ylabel('$f_{in}$ (Hz)', 'Interpreter', 'latex')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca, 'TickLength', [0.007, 0.007]) 
xticks(datenum(t_ticks_dm))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([DM.time(2), DM.time(end)])
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca,'Box','on')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu', 100)))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{in}^{5} \cdot$ S($t,f_{in}$) $\cdot \frac{df_{ob}}{df_{in}}$ (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [ 10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ;  
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

%------------ Subplot 4 ------------% 
subplot(2,2,4)

% Plot intrinsic frequency saturation spectrogram from SMODE2021
pc = pcolor(SM.time, SM.f(Ifreq), SM.spec_in(Ifreq,:));

% Set figure attributes
title('(d)')
pc.EdgeColor = 'none';
xlim([SM.time(1), SM.time(end)])
xlabel('UTC time from Oct 29$^{\textrm{th}}$, 2021', 'Interpreter', 'latex')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca, 'TickLength', [0.007, 0.007]) 
xticks(datenum(t_ticks_sm))
datetick('x', 'mmm dd', 'keepticks')
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca,'Box','on')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu', 100)))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{in}^{5} \cdot$ S($t,f_{in}$) $\cdot \frac{df_{ob}}{df_{in}}$ (m$^2$ Hz$^{4}$)';
caxis([10^-7, 2*10^-3]);
cb.Ticks = [10^-7; 10^-6; 10^-5; 10^-4; 10^-3];
cb.TickLabels = {'$10^{-7}$'; '$10^{-6}$'; '$10^{-5}$'; '$10^{-4}$'; '$10^{-3}$'};  
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

% Save Figure
saveas(gcf, [fig_path 'figure_10.png'])
