function plot_dir_spec(freq, theta, Sd, Nspokes, Ncircles, Rticks, Contours, cb_ticks, cb_ticklabels, fontsize, cb_label)

    %%%%
    % plot_dir_spec(freq, theta, Sd, Nspokes, Ncircles, Rticks, Contours, cb_ticks, cb_ticklabels, fontsize, cb_label)
    %
    % Function for plotting a contourf directional wave spectrum in
    % frequency and wave direction space (frequency axis is logarithmic).  
    %
    %   Parameters
    %   ---------- 
    %   freq     : Wave cyclical frequency array. 
    %   theta    : Wave direction array. 
    %   Sd       : Directional wave spectrum. 
    %   Nspokes  : Number of radial lines corresponding to the theta axis
    %              tick marks. 
    %   Ncircles : Number of azimuthal circles corresponding to the
    %              freqeuncy axis tick marks
    %   Rticks   : Cell array of frequency labels. I usually obtain these
    %              values by obtained by the following: 
    %               logspace(log10(freq(2)), log10(freq(end)), Ncircles))  
    %   Contours : 
    %   cb_ticks : Array of doubles specifying the logarithmic colorbar 
    %              tick mark locations for the power spectral density.  
    %   cb_ticklabels : Cell array of strings specifying the logarithmic
    %                   colorbar ticklabels. Must be the same size as the
    %                   cb_ticks array.
    %   cb_label : String specifying the colorbar's label. 
    %   fontsize : Fontsize of text. 
    % 
    %   Returns
    %   -------
    %   Directional wave spectrum figure. 
    %   
    %%%%
    
    % Set plotting parameters
    pos = linspace(freq(2), freq(end), Ncircles);                           % Position of frequency labels
    theta_s = linspace(0,360,length(theta));                                % Set theta for plotting                                                  
    cmap = colormap(flipud(cbrewer2('RdYlBu', numel(Contours)))); 
    Contours_n = [10^-16, Contours];
    cmap_n = cat(1,cmap(1,:), cmap);
    LineWidth = 0.5;
    LineColor = 'k';
    LineStyle = '-';
    TextColor = 'w';

    % Plot Directional spectra  
    [~,cb] = polarcontourf(freq(2:end), theta_s, Sd(2:end,:), ...
                           'circlesPos', pos, 'Rscale', 'log',...
                           'Ncircles', Ncircles, 'Nspokes', Nspokes,...
                           'typeRose' ,'meteo', 'Contours', Contours_n,...
                           'colBar', 3, 'fontsize', fontsize, 'colormap', cmap_n,...
                           'LineWidth', LineWidth, 'LineColor', LineColor,...
                           'LineStyle', LineStyle, 'RtickLabel', Rticks, ...
                           'TextColor', TextColor);

    % Set logarithmic colorbar
    cb.Ticks = log10(cb_ticks);
    cb.TickLabels = cb_ticklabels; 
    cb.Label.String = cb_label; 'S$_{ob}$($f_{ob}$, $\theta$) (m$^2$ Hz$^{-1}$ deg$^{-1}$)';
    cb.TickDirection = 'out';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fontsize;
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontsize;
    
    % Set other figure attributes 
    set(gcf,'color','w')
    set(gca,'FontSize',fontsize)
    set(gca,'TickLabelInterpreter','latex')
    
end