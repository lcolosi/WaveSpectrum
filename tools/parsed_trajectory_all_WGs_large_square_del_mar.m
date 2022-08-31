%% Parse Trajectories of wave gliders for large square during delmar deployment
% Luke Colosi | lcolosi@ucsd.edu | April 5th, 2021
% Code adapted from code created by Laurent Grare 

clc, clear, close all

%%%%%%%%% Notes %%%%%%%%%
% (1) String Angles 000: Right leg, 090: bottom leg, 180: left leg, 270: top leg
% (2) Kelvin Large box legs are not great because the PLD2_1hz_ALL.mat file
% is missing a lot of data (does not include any data from the small box)
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Initial global variables %%%%%%%%%%
runmean_time = 60;                                                          % window length for running mean 
DELTA = 25;                                                                 % Offset distance from the vertices of WG trajectory 
LEG = {'000';'090';'180';'270'};                                            % String for wave glider legs

% Set path to data
%%%% Vehicle %%%%
% vehicle = 'PLANCK';
vehicle = 'STOKES';
% vehicle = 'KELVIN';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';

% Load data
load([ROOT vehicle '/PLD2_1hz_ALL.mat'])
SLOW  = PLD2_1hz;

%% Process data

% Set mean positions and vertices of rectangle#1 (Large rectangle)
if strcmp(vehicle, 'PLANCK')
    LON0 = -117.675;  LAT0 = 32.995;  % Mean Position
    LON1 = -117.67;   LAT1 = 33.00;  % Vertex 1
    LON2 = -117.68;   LAT2 = 33.00;  % Vertex 2
    LON3 = -117.68;   LAT3 = 32.99;  % Vertex 3
    LON4 = -117.67;   LAT4 = 32.99;  % Vertex 4
elseif strcmp(vehicle, 'STOKES')
    LON0 = -117.645;  LAT0 = 32.995;
    LON1 = -117.64;   LAT1 = 33.00;
    LON2 = -117.65;   LAT2 = 33.00;
    LON3 = -117.65;   LAT3 = 32.99;
    LON4 = -117.64;   LAT4 = 32.99;
elseif strcmp(vehicle, 'KELVIN')
    LON0 = -117.66;  LAT0 = 32.995;
    LON1 = -117.655;   LAT1 = 33.00;
    LON2 = -117.665;   LAT2 = 33.00;
    LON3 = -117.665;   LAT3 = 32.99;
    LON4 = -117.655;   LAT4 = 32.99;
end

% Compute distance between navigational waypoints for top (270 degrees) and left (180 degrees) legs
[~,Lx] = legs([LAT1 LAT2],[LON1 LON2]);
[~,Ly] = legs([LAT2 LAT3],[LON2 LON3]);

% Convert to units from nautical miles to meters 1 Nautical mile = 1852
% meters
Lx = Lx * 1852;
Ly = Ly * 1852;

% Set offset distance from the center of WG trajectory
Dx = 75;
Dy = 75;

% Compute latitude at minutes resolution from degrees (integer) and minutes
% (1 degree = 60 minutes)
SLOW.latitude  = SLOW.latitude_d_1hz  + SLOW.latitude_m_1hz/60;
SLOW.longitude = SLOW.longitude_d_1hz + SLOW.longitude_m_1hz/60;

% Compute running mean of longitude and latitude
SLOW.runmean.longitude = movmean(SLOW.longitude,60,'omitnan');
SLOW.runmean.latitude  = movmean(SLOW.latitude,60,'omitnan');

% Convert running mean longitude and latitude from degrees to meters with
% (meters corresponds to distance from center of box)
[SLOW.runmean.Xe,SLOW.runmean.Yn] = lonlat2m(SLOW.runmean.longitude,SLOW.runmean.latitude,LON0,LAT0);

% Compute heading and distance between each of the running mean longitude 
% and latitude waypoints
[SLOW.runmean.course,SLOW.runmean.dist] = legs(SLOW.runmean.latitude,SLOW.runmean.longitude);

% Find indices of data points for each leg 
I000 = find(SLOW.runmean.Xe < (+Lx+Dx)/2 & SLOW.runmean.Xe > (+Lx-Dx)/2 & SLOW.runmean.Yn > (-Ly-Dy)/2 & SLOW.runmean.Yn < (+Ly+Dy)/2);  % 
I090 = find(SLOW.runmean.Xe > (-Lx-Dx)/2 & SLOW.runmean.Xe < (+Lx+Dx)/2 & SLOW.runmean.Yn > (-Ly-Dy)/2 & SLOW.runmean.Yn < (-Ly+Dy)/2);
I180 = find(SLOW.runmean.Xe > (-Lx-Dx)/2 & SLOW.runmean.Xe < (-Lx+Dx)/2 & SLOW.runmean.Yn > (-Ly-Dy)/2 & SLOW.runmean.Yn < (+Ly+Dy)/2);
I270 = find(SLOW.runmean.Xe > (-Lx-Dx)/2 & SLOW.runmean.Xe < (+Lx+Dx)/2 & SLOW.runmean.Yn > (+Ly-Dy)/2 & SLOW.runmean.Yn < (+Ly+Dy)/2);

% Remove indices near the vertices of legs
I000 = intersect(I000,find_wrap(SLOW.runmean.course,000-DELTA ,000+DELTA,360) );
I090 = intersect(I090,find_wrap(SLOW.runmean.course,090-DELTA ,090+DELTA,360) );
I180 = intersect(I180,find_wrap(SLOW.runmean.course,180-DELTA ,180+DELTA,360) );
I270 = intersect(I270,find_wrap(SLOW.runmean.course,270-DELTA ,270+DELTA,360) ); 

% Compile all leg indicies 
Ileg = union(I000,union(I090,union(I180,I270)));

% Create a vector the length of observations 
Iall = 1:length(SLOW.time);

% Find indices for data points not categorized as a leg
Irest = setdiff(Iall,Ileg);

%% Index data from Legs and Save data

% Loop through legs
for j = 1:4
    
    % Call string for leg
    leg = LEG{j};
    
    % Initialize indice for legs
    eval(['I = I' leg ';'])
    
    % Find end indices of leg segments
    Je = find(diff(I)/60 > 20);
   
    % Initialize start and end indices of leg segments vectors
    Js = Je + 1; 
    Je = [Je ; length(I)];
    Js = [1 ; Js];
    
    % Plot 
    clf, close all
    %%%% Subplot 1 %%%%
    subplot(121)
        hold on; grid on;
        plot(I,'.')
        plot(Js,I(Js),'or')
        plot(Je,I(Je),'og')
    %%%% Subplot 2 %%%%
    subplot(122)
        hold on; grid on;
        plot(SLOW.runmean.Xe,SLOW.runmean.Yn,'.b');
        axis square
        axis equal
        xlim([-Lx/2-Dx Lx/2+Dx])
        ylim([-Ly/2-Dy Ly/2+Dy])
        
    % Set fieldnames to loop through and counter
    field = fieldnames(SLOW);
    m = 0;
    
    % Loop through starting indices
    for i = 1:length(Js)
        
        % Find time indices for ith leg
        It = find(SLOW.time >= SLOW.time(I(Js(i))) & SLOW.time <= SLOW.time(I(Je(i))));
        
        % Compute distance traversed by ith leg
        [~,length_leg] = legs(SLOW.latitude(I([Js(i) Je(i)])) , SLOW.longitude(I([Js(i) Je(i)])));  % units: nautical miles
        
        % Convert from nautical miles to meters
        length_leg = length_leg * 1852;
        
        % Compute the length of data points in ith leg
        ni = length(It);
        
        % Discard leg if number of data points is less than 1200 or
        % distance traversed during leg is less than 250 meters
        if ni < 1200 | length_leg < 250
            continue
        else
            % Reinitialize counter
            m = m + 1;
            
            % Save starting and ending indices and times for ith leg
            Ks(m) = I(Js(i));
            Ke(m) = I(Je(i));
            Ts(m) = SLOW.time(It(1));
            Te(m) = SLOW.time(It(end));
        end
        
        % Plot ith leg
        %%%% Subplot 1 %%%%
        subplot(121)
            plot(Js(i):Je(i),I(Js(i):Je(i)),'xk')
        %%%% Subplot 2 %%%%
        subplot(122)
            plot(SLOW.runmean.Xe(It),SLOW.runmean.Yn(It),'.r');
    end
    
    % Save Start and End time of leg 
    eval(['LEG_' leg '.time_start_leg = Ts;'])
    eval(['LEG_' leg '.time_end_leg = Te;'])
    
end

% Save start and end time of legs in mat file
save([ROOT 'LEGS_' vehicle '_LARGE_SQUARE.mat'],'LEG_000','LEG_090','LEG_180','LEG_270')

%% Plot Legs

clf, close all

% Set colorbar
coul = jet(4);

% Plot trajectories
subplot(1,2,1)

    hold on; 
    %plot(SLOW.runmean.Xe,SLOW.runmean.Yn,'.m')
    plot(SLOW.runmean.Xe(I000),SLOW.runmean.Yn(I000),'.','color',coul(1,:))
    plot(SLOW.runmean.Xe(I090),SLOW.runmean.Yn(I090),'.','color',coul(2,:))
    plot(SLOW.runmean.Xe(I180),SLOW.runmean.Yn(I180),'.','color',coul(3,:))
    plot(SLOW.runmean.Xe(I270),SLOW.runmean.Yn(I270),'.','color',coul(4,:))
    hold off

    % Set figure Attributes
    axis square
    axis equal
    xlim([-Lx/2-Dx Lx/2+Dx])
    ylim([-Ly/2-Dy Ly/2+Dy])
    grid on

% Smooth COG_1hz with low pass filter 
SLOW.true_course = mod(atan2d(movmean(sind(180 + SLOW.COG_1hz), 60*4), movmean(cosd(180 + SLOW.COG_1hz), 60*4)),360); %convention: CW, going towards, ref north 
    
% Plot True Course
subplot(1,2,2)

    plot(SLOW.time(I000), SLOW.true_course(I000),'.','color',coul(1,:))
    hold on
    plot(SLOW.time(I090), SLOW.true_course(I090),'.','color',coul(2,:))
    plot(SLOW.time(I180), SLOW.true_course(I180),'.','color',coul(3,:))
    plot(SLOW.time(I270), SLOW.true_course(I270),'.','color',coul(4,:))
    hold off
    
    % Set figure Attributes
    axis square
    datetick('x', 'mm/dd')
    grid on