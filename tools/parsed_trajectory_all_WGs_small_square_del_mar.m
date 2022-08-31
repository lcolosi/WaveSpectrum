%% Parse Trajectories of wave gliders for small square during delmar deployment
% Luke Colosi | lcolosi@ucsd.edu | April 5th, 2021
% Code adapted from Code created by Laurent Grare 

clc, clear, close all

%%%%%%%%% Notes %%%%%%%%%
% (1) 000: Right leg, 090: bottom leg, 180: left leg, 270: top leg
% (2) Kelvin Large box legs are terrible because the PLD2_1hz_ALL.mat file
% is missing a lot of data (does not include any data from the small box)
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Initial global variables %%%%%%%%%%
runmean_time = 60;                                                          % wind length for running mean 
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
if strcmp(vehicle, 'KELVIN')
    LON1 = -117.565;  LAT1 = 32.91;  % Vertex 1
    LON2 = -117.565;  LAT2 = 32.905;  % Vertex 2
    LON3 = -117.57;  LAT3 = 32.905;  % Vertex 3
    LON4 = -117.57;  LAT4 = 32.91;  % Vertex 4
    LON0 = -117.5675;  LAT0 = 32.9075;  % Mean Position
elseif strcmp(vehicle, 'STOKES')
    LON1 = -117.56075;  LAT1 = 32.91750;
    LON2 = -117.56075;  LAT2 = 32.91250;
    LON3 = -117.56575;  LAT3 = 32.91250;
    LON4 = -117.56575;  LAT4 = 32.91750;
    LON0 = -117.56325;  LAT0 = 32.915;
elseif strcmp(vehicle, 'PLANCK')
    LON1 = -117.56925;  LAT1 = 32.91750;
    LON2 = -117.56925;  LAT2 = 32.91250;
    LON3 = -117.57425;  LAT3 = 32.91250;
    LON4 = -117.57425;  LAT4 = 32.91750;
    LON0 = -117.57175;  LAT0 = 32.915;
end

% Compute distance between navigational waypoints for top (270 degrees) and left (180 degrees) legs
[~,Lx] = legs([LAT2 LAT3],[LON2 LON3]);
[~,Ly] = legs([LAT1 LAT2],[LON1 LON2]);

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
I000 = find(SLOW.runmean.Xe > (-Lx-Dx)/2 & SLOW.runmean.Xe < (-Lx+Dx)/2 & SLOW.runmean.Yn > (-Ly-Dy)/2 & SLOW.runmean.Yn < (+Ly+Dy)/2);
I090 = find(SLOW.runmean.Xe > (-Lx-Dx)/2 & SLOW.runmean.Xe < (+Lx+Dx)/2 & SLOW.runmean.Yn > (+Ly-Dy)/2 & SLOW.runmean.Yn < (+Ly+Dy)/2);
I180 = find(SLOW.runmean.Xe > (+Lx-Dx)/2 & SLOW.runmean.Xe < (+Lx+Dx)/2 & SLOW.runmean.Yn > (-Ly-Dy)/2 & SLOW.runmean.Yn < (+Ly+Dy)/2);
I270 = find(SLOW.runmean.Xe > (-Lx-Dx)/2 & SLOW.runmean.Xe < (+Lx+Dx)/2 & SLOW.runmean.Yn > (-Ly-Dy)/2 & SLOW.runmean.Yn < (-Ly+Dy)/2);

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
        if ni < 600 | length_leg < 250
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
save([ROOT 'LEGS_' vehicle '_SMALL_SQUARE.mat'],'LEG_000','LEG_090','LEG_180','LEG_270')

%% Plot Legs

clf, close all

% Set colorbar
coul = jet(4);

% Plot trajectories
hold on; 
plot(SLOW.runmean.Xe,SLOW.runmean.Yn,'.m')
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
