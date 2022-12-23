%% Figure 3: Mapping observed to intrinsic frequency
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2021

%-------------------------------- Caption --------------------------------%
% Mapping from observed to intrinsic frequency for a platform moving at 
% 0.5 ms$^{-1}$ (a), 1 ms$^{-1}$ (b), and 2 ms$^{-1}$ (c) with line color
% corresponding to the platform's direction relative to the waves. The
% platform is moving with the waves when $\theta_r$ ranges from 
% $0^\circ \leq \theta_r < 90^\circ$ (red curves) and is moving against 
% the waves when $\theta_r$ ranges from 
% $90^\circ < \theta_r \leq 180^\circ$ (blue curves). The platform is
% moving perpendicular to the waves when $\theta_r = 90^\circ$. The gray
% dashed line is the cutoff frequency $f_c$ with gray circles being $f_c$
% for a given directions. The black dashed line is the one-to-one line
% occurring when the projected speed vanishes. 
%-------------------------------------------------------------------------%

clc, clear, close all;

% Set text interpreter 
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

% Set path for figures
fig_path = '../figs/'; 

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
g = 9.81;                                                                   % Gravitational acceleration (units: ms^-2)
U = [0.5, 1, 2];                                                            % Speed of wave glider (units: ms^-1)
theta_r = (0:30:180);                                                       % Relative angle between mean platform propagation and wave direction (units: degrees, reference north, CW, going towards)

% Frequency parameters
df_ob = 0.001;                                                              % Observed frequency resolution (units: Hz)
df_in = 0.001;                                                              % Intrinsic frequency resolution (units: Hz)
f_ob = (0:df_ob:1)';                                                        % Cyclical observed frequency (units: Hz)
f_in = (0:df_in:1)';                                                        % Cyclical intrinsic frequency (for inverse mapping f_ob(f_in)) (units: Hz)
f_cut = 1;                                                                  % Noise cutoff frequency (units: Hz)
F.f_in = nan(length(f_ob), length(U), length(theta_r));                     % Mapped cyclical observed frequency 
F.f_ob = nan(length(f_ob), length(U), length(theta_r));                     % Mapped cyclical intrinsic frequency 
F.deriv_in = nan(length(f_ob), length(U), length(theta_r));                 % Mapped first derivative of intrinsic frequency
F.f_in_max = nan(length(U), length(theta_r));                               % Intrinsic frequency from fowards mapping f_in(f_ob) where df_in/df_ob approaches infinity
F.f_ob_max = nan(length(U), length(theta_r));                               % Observed frequency from fowards mapping f_in(f_ob) where df_in/df_ob approaches infinity

%% Compute intrinsic cyclical frequency

% Loop through speeds
for ispeed = 1:length(U)
    
    % Call speed
    iU = U(ispeed);
    
    % Loop through relative angle
    for iangle = 1:length(theta_r)
        
        % Call relative angle
        itheta = theta_r(iangle);
        
        %%%%%%%%%% Branch Conditional statements %%%%%%%%%%
        %---------- Branch 1 ----------%
        if cosd(itheta) < 0 && iU > 0
            disp(['U = ', num2str(iU), ', theta_r = ', num2str(itheta), ', Branch 1'])
            
            % Compute intrinsic frequency for branch 1
            F.f_in(1:length(f_ob),ispeed,iangle) = f_ob + (g - 4*pi*f_ob*iU*cosd(itheta) - sqrt(g^2 - 8*pi*f_ob*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
            
            % Compute first derivative of intrinsic frequency
            F.deriv_in(1:length(f_ob),ispeed,iangle) = -((g)/(sqrt(g^2 - 8*pi*g*iU*cosd(itheta)*f_ob)));
            
            % Save observed frequency
            F.f_ob(:,ispeed,iangle) = f_ob;
            
        %---------- Branch 2 ----------%
        elseif cosd(itheta) == 0 || iU == 0
            disp(['U = ', num2str(iU), ', theta_r = ',num2str(itheta), ', Branch 2'])
            
            % Compute intrinsic frequency for branch 2 
            F.f_in(1:length(f_ob),ispeed,iangle) = f_ob; 
            
            % Compute first derivative of intrinsic frequency
            F.deriv_in(1:length(f_ob),ispeed,iangle) = df_ob * ones(length(f_ob), 1);
            
            % Save observed frequency
            F.f_ob(:,ispeed,iangle) = f_ob;
            
        %---------- Branch 3-5 ----------%
        elseif cosd(itheta) > 0 && iU > 0
            disp(['U = ', num2str(iU), ', theta_r = ', num2str(itheta), ', Branch 3-5'])
            
            % Compute f_ob where df_in(f_ob)/df_ob tends towards infinity  
            F.f_ob_max(ispeed,iangle) = (g)/(8*pi*iU*cosd(itheta));
            
            % Compute f_in at the f_ob value where df_in(f_ob)/df_ob tends towards infinity
            F.f_in_max(ispeed,iangle) = (g)/(4*pi*iU*cosd(itheta));
            
            % Compute f_in when f_ob = 0
            F.f_in_zero(ispeed,iangle) = (g)/(2*pi*iU*cosd(itheta));
            
            %---------- Platform moving slower than energy and crests ----------%
            if F.f_in_max(ispeed,iangle) > 1
                disp('Branch 3: Platform moving slower than energy and crests')
                
                % Compute the f_ob where f_in = f_cut for branch 3
                f_ob_cut = -((2*pi*iU*cosd(itheta)*f_cut^2)/(g)) + f_cut;
                f_in_cut = f_ob_cut + (g - 4*pi*f_ob_cut*iU*cosd(itheta) - sqrt(g^2 - 8*pi*f_ob_cut*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                
                % Obtain f_ob values that map to f_in less than or equal to f_ob_cut
                f_ob_map = f_ob(f_ob <= f_ob_cut);
                
                % Compute intrinsic frequency for branch 3
                F.f_in(1:length(f_ob_map),ispeed,iangle) = f_ob_map + (g - 4*pi*f_ob_map*iU*cosd(itheta) - sqrt(g^2 - 8*pi*f_ob_map*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                
                % Compute first derivative of intrinsic frequency
                F.deriv_in(1:length(f_ob_map),ispeed,iangle) = -((g)/(sqrt(g^2 - 8*pi*g*iU*cosd(itheta)*f_ob_map)));
                
                % Save observed frequency
                F.f_ob(1:length(f_ob_map),ispeed,iangle) = f_ob_map;
            
            %---------- Platform moving faster than energy but slower than crests ----------%
            elseif F.f_in_max(ispeed,iangle)  > 0 && F.f_in_max(ispeed,iangle) < 1 && F.f_in_zero(ispeed,iangle) > 1
                disp('Branch 3 and 4: Platform moving faster than energy but slower than crests')
                
                % Compute the f_ob where f_in = f_cut for branch 4
                f_ob_cut = -((2*pi*iU*cosd(itheta)*f_cut^2)/(g)) + f_cut;
                f_in_cut = f_ob_cut + (g - 4*pi*f_ob_cut*iU*cosd(itheta) + sqrt(g^2 - 8*pi*f_ob_cut*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                
                % Obtain f_ob values for mapping to branch 3 (f_ob less than or equal to f_ob_max)
                f_ob_map3 = f_ob(f_ob <= F.f_ob_max(ispeed,iangle));
                
                % Obtain f_ob values for mapping to branch 4
                f_ob_map4 = f_ob_map3(f_ob_map3 >= f_ob_cut);
                
                % Compute intrinsic frequency for branch 3 and 4
                f_in_b3 = f_ob_map3 + (g - 4*pi*f_ob_map3*iU*cosd(itheta) - sqrt(g^2 - 8*pi*f_ob_map3*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                f_in_b4 = f_ob_map4 + (g - 4*pi*f_ob_map4*iU*cosd(itheta) + sqrt(g^2 - 8*pi*f_ob_map4*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                
                % Compute first derivative of intrinsic frequency for branches 3 and 4
                deriv_in_3 = -((g)/(sqrt(g^2 - 8*pi*g*iU*cosd(itheta)*f_ob_map3)));
                deriv_in_4 = ((g)/(sqrt(g^2 - 8*pi*g*iU*cosd(itheta)*f_ob_map4)));
                
                % Combine f_in from branches 3 and 4
                F.f_in(1:length([f_ob_map3; f_ob_map4]),ispeed,iangle) = [f_in_b3; flip(f_in_b4)];
                
                % Combine deriv_f_in from branches 3 and 4
                F.deriv_in(1:length([f_ob_map3; f_ob_map4]),ispeed,iangle) = [deriv_in_3, deriv_in_4];
                
                % Save observed frequency
                F.f_ob(1:length([f_ob_map3; f_ob_map4]),ispeed,iangle) = [f_ob_map3; flip(f_ob_map4)];
                
            %---------- Platform moving faster than energy and crests ----------%
            elseif F.f_in_zero(ispeed,iangle) < 1
                disp('Branch 3, 4, and 5 : Platform moving faster than energy and crests')
                
                % Compute the f_ob where f_in = f_cut for branch 5
                f_ob_cut = ((2*pi*iU*cosd(itheta)*f_cut^2)/(g)) - f_cut;
                f_in_cut = -f_ob_cut + (g + 4*pi*f_ob_cut*iU*cosd(itheta) + sqrt(g^2 + 8*pi*f_ob_cut*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                
                % Obtain f_ob values for mapping to branch 3 and 4 (f_ob less than or equal to f_ob_max)
                f_ob_map_3_4 = f_ob(f_ob <= F.f_ob_max(ispeed,iangle));
                disp([num2str(F.f_ob_max(ispeed,iangle)) ', ' num2str(f_ob_map_3_4(end))])
                
                % Obtain f_ob values for mapping to branch 5
                f_ob_map5 = f_ob(f_ob <= f_ob_cut);
                
                % Compute intrinsic frequency for branch 3, 4, and 5
                f_in_b3 = f_ob_map_3_4 + (g - 4*pi*f_ob_map_3_4*iU*cosd(itheta) - sqrt(g^2 - 8*pi*f_ob_map_3_4*g*iU*cosd(itheta)))/(4*pi*iU * cosd(itheta));
                f_in_b4 = f_ob_map_3_4 + (g - 4*pi*f_ob_map_3_4*iU*cosd(itheta) + sqrt(g^2 - 8*pi*f_ob_map_3_4*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                f_in_b5 = -f_ob_map5 + (g + 4*pi*f_ob_map5*iU*cosd(itheta) + sqrt(g^2 + 8*pi*f_ob_map5*g*iU*cosd(itheta)))/(4*pi* iU * cosd(itheta));
                
                % Compute first derivative of intrinsic frequency for branches 3 and 4
                deriv_in_3 = -((g)/(sqrt(g^2 - 8*pi*g*iU*cosd(itheta)*f_ob_map_3_4)));
                deriv_in_4 = ((g)/(sqrt(g^2 - 8*pi*g*iU*cosd(itheta)*f_ob_map_3_4)));
                deriv_in_5 = ((g)/(sqrt(g^2 + 8*pi*g*iU*cosd(itheta)*f_ob_map5)));
                
                % Combine branches  
                F.f_in(1:length([f_in_b3; f_in_b4; f_in_b5]),ispeed,iangle) = [f_in_b3; flip(f_in_b4); f_in_b5];
                
                % Combine deriv_f_in from branches 3 and 4
                F.deriv_in(1:length([f_ob_map_3_4; f_ob_map_3_4; f_ob_map5]),ispeed,iangle) = [deriv_in_3, deriv_in_4, deriv_in_5];
            
                % Save observed frequency
                F.f_ob(1:length([f_ob_map_3_4; f_ob_map_3_4; f_ob_map5]),ispeed,iangle) = [f_ob_map_3_4; flip(f_ob_map_3_4); -f_ob_map5];
                
            end
        end
        
        % Compute the inverse mapping (f_in to f_ob)
        F.f_ob_im(:,ispeed,iangle) = f_in - ((2*pi*iU*cosd(itheta))/(g)) * f_in.^2;
        
        % Compute the first derivative of f_ob wrt f_in
        F.deriv_ob_im(:,ispeed,iangle) = 1 - ((4*pi*iU*cosd(itheta))/(g)) * f_in;
        
        % Compute the frequency resolution of f_in
        F.df_ob_im(:,ispeed,iangle) = F.deriv_ob_im(:,ispeed,iangle)*df_in;
        
    end
end

%% Plot Observed vs Intrinsic Frequency 
clc, close all 

% Create Figure and axes
fig = figure( 'Name', 'Observed vs Intrinsic Frequency');

% Set figure attributes
POS = [100 100 500 2000];                                                   %[100 100 2000 300];
set(gcf,'color',[1 1 1])
set(gcf,'position',POS) 
fontsize = 13;

% Obtain RGB triplet for colormap
cmap = colormap(cbrewer2('RdYlBu'));

% Bin average colorbar to match number of relative angle cases
ntriplet = floor(length(cmap)/length(theta_r))*length(theta_r);             % number of RGB triplets that divide evenly into number of angle bins
ncolor = floor(length(cmap)/length(theta_r));                               % number of triplet to be bin averaged
cmap_bin = sepblockfun(cmap(1:ntriplet,:), [ncolor,1], @mean);

% Set colormap and colorbar variables
z_m = linspace(0,180,7);
z = 0:12.8571:180;

% Obtain RGB triplet of colormap
cmap = colormap(cmap_bin);

% Plot intrinsic frequency as a function of observed frequency
for isubplot = 1:length(U)
    
    % Initialize subplot
    subplot(3,1,isubplot)
    
    % Loop through angles
    for iangle = 1:length(theta_r)
        
        % Set line color
        [~,idx_c] = min(abs(z_m-theta_r(iangle)));
        
        hold on 
        
            %Plot f_in as a function of f_ob in linear space
            plot(F.f_ob(:,isubplot,iangle), F.f_in(:,isubplot,iangle), '-', 'LineWidth', 2, 'Color', cmap(idx_c,:))

            % Plot Ucos(theta_r) = 0 and Ucos(theta_r) = -c_g 
            pc1 = plot(f_ob, f_ob, '--k', 'LineWidth', 2);
            pc2 = plot(f_ob, (2)*f_ob, 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

            % Plot f_in at bifuraction frequencies
            if F.f_in_max(isubplot,iangle) < 1 && F.f_in_max(isubplot,iangle) > 0
                pc3 = plot(F.f_ob_max(isubplot,iangle), F.f_in_max(isubplot,iangle), 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 12, 'LineWidth', 2);
            end

        hold off

    end

    % Set figure attributes
    axis equal 
    axis tight
    if isubplot == 1
        ax1 = gca; 
        label = '(a)';
    elseif isubplot == 2  
        label = '(b)'; 
    elseif isubplot == 3
        ax3 = gca; 
        label = '(c)';
        xlabel('$f_{ob}$ (Hz)','Interpreter','Latex')
        leg = legend([pc1, pc2, pc3], 'U$\cos(\theta_r) = 0$', 'U$\cos(\theta_r) = -c_g$', '$f_c(f_{ob})$','location', 'northwest', 'fontsize', fontsize);
    end
    title([label ' U $= \;$', num2str(U(isubplot)), ' ms$^{-1}$'])
    ylabel('$f_{in}$ (Hz)','Interpreter','Latex')
    xlim([-0.4,1])
    ylim([0,1])
    xticks(-0.4:0.2:1)
    yticks(0:0.2:1)
    grid on
    set(gca,'FontSize',fontsize)
    set(gca,'Box','on')

end

%-------- Set Colorbar --------%
colormap(cmap_bin);
cb = colorbar;
cb.Position = [0.8201 0.3949 0.0279 0.2415]; 
cb.Ticks = z(2:2:end)/max(z);
cb.TickLabelInterpreter = 'Latex'; 
cb.Label.Interpreter = 'Latex';
cb.TickLabels = {'0', '30', '60', '90', '120', '150', '180'};   
cb.FontSize = fontsize;
cb.TickLabels = num2str(theta_r');
cb.Label.String = '';
cb.TickDirection = 'out';
annotation('textbox',[0.8130 0.6348 0.1080 0.0359],...
           'String',{'$\theta_r \;(^\circ)$'},'Interpreter','latex',...
           'FontSize',fontsize,'EdgeColor','none');

%-------- Set poition of legend --------%
set(leg,'Position',[0.303 0.87889 0.2689 0.0042])

%-------- Set the position of the subplots --------%
% find current position [x,y,width,height]
pos1 = get(ax1,'Position');
pos3 = get(ax3,'Position');

% set width and height of third axis equal to first
pos3(3) = pos1(3);
pos3(4) = pos1(4);
set(ax3,'Position',pos3)

% Save Figure
saveas(fig, strcat(fig_path, 'figure_3.png'))
