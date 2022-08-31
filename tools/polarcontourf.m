function [varargout] = polarcontourf(R,theta,Z,varargin)
% [h,c] = polarcontourf(R,theta,Z,varargin) is a contour plot of matrix 
% Z for a vector radius R and a vector angle theta. 
% The elements of Z specify the color in each cell of the 
% plot. The goal is to apply contourf function with a polar grid, which 
% provides a better visualization than a cartesian grid.
%
%% Syntax
% 
% [h,c] = polarcontourf(R,theta,Z,)
% [h,c] = polarcontourf(R,theta,Z,'Ncircles',10)
% [h,c] = polarcontourf(R,theta,Z,'Nspokes',5)
% [h,c] = polarcontourf(R,theta,Z,'Nspokes',5,'colBar',0) 
% [h,c] = polarcontourf(R,theta,Z,'Nspokes',5,'labelR','r (km)')
% 
% INPUT
%	* R :
%        - type: float
%        - size: [1 x Nrr ] where Nrr = numel(R).
%        - dimension: radial distance.
%	* theta : 
%        - type: float
%        - size: [1 x Ntheta ] where Ntheta = numel(theta).
%        - dimension: azimuth or elevation angle (deg).
%        - N.B.: The zero is defined with respect to the North.
%	* Z : 
%        - type: float
%        - size: [Ntheta x Nrr]
%        - dimension: user's defined .
%	* varargin:
%        - Ncircles: number  of circles for the grid definition.
%        - autoOrigin: 'on' (the first circle of the plar grid has a radius
%          equal to the lowest value of R) or 'off'.
%        - Nspokes: number of spokes for the grid definition.
%        - colBar: display the colorbar or not. Options include: 
%                  colBar = 1: Linear colorbar
%                  colBar = 2: Logarithmic colorbar (no lower or upper 
%                              color set for values outside levels range)
%                  colBar = 3: Logarithmic colorbar with lower color set
%                              for values outside levels range. 
%                  Any other value other than 1, 2, or 3 does not generate
%                  a colorbar. 
%        - labelR: title for radial axis.
%        - RtickLabel: Tick label for the radial axis.
%        - colormap: Colormap for the contourf function
%        - ncolor: Number of colors in the colorbar and contourf
%        - circlesPos: position of the circles with respect to the origin (it overwrites Ncircles if necessary)
%        - LineStyle: Line style for axes. 
%        - LineWidth: Line width for axes. 
%        - Line Color: Color of axes. 
% 
% 
% OUTPUT
% h: returns a handle to a SURFACE object.
% c: returns a handle to a COLORBAR object.
%
%% Examples 
% R = linspace(3,10,100);
% theta = linspace(0,180,360);
% Z = linspace(0,10,360)'*linspace(0,10,100);
%
% figure
% polarcontourf(R,theta,Z,'Contours', Contours, 'Ncircles',3)
%
%% Author
% Etienne Cheynet, University of Stavanger, Norway. 23/10/2019 and Luke
% Colosi, University of California, San Diego, USA. 15/08/2022 
% see also contourf
% 

%%  InputParseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Ncircles',5);
p.addOptional('autoOrigin','on');
p.addOptional('Nspokes',8);
p.addOptional('labelR','');
p.addOptional('RtickLabel',[]);
p.addOptional('colBar',1);
p.addOptional('Contours',[]);
p.addOptional('Rscale','linear');
p.addOptional('colormap','parula');
p.addOptional('ncolor',[]);
p.addOptional('typeRose','default'); % 'meteo' or 'default'
p.addOptional('circlesPos',[]);
p.addOptional('fontsize',8);
p.addOptional('LineStyle', '-')
p.addOptional('LineColor', 'k')
p.addOptional('LineWidth', 1)
p.addOptional('TextColor', 'k')
p.parse(varargin{:});

Ncircles = p.Results.Ncircles;
Nspokes = p.Results.Nspokes ;
labelR = p.Results.labelR ;
RtickLabel = p.Results.RtickLabel ;
colBar = p.Results.colBar ;
Contours = p.Results.Contours;
Rscale = p.Results.Rscale ;
autoOrigin = p.Results.autoOrigin ;
myColorMap = p.Results.colormap ;
ncolor = p.Results.ncolor ;
circPos = p.Results.circlesPos ;
typeRose = p.Results.typeRose ;
fontsize = p.Results.fontsize ;
LineStyle = p.Results.LineStyle;
LineColor = p.Results.LineColor;
LineWidth = p.Results.LineWidth;
TextColor = p.Results.TextColor;

if isempty(Contours)
    error(' Contours must be defined!!! ')
end

if ~isempty(circPos)
    Origin = max([min(circPos),min(R)]);
    circPos(circPos<min(R))=[];
    circPos(circPos>max(R))=[];
elseif strcmpi(autoOrigin,'on')
    Origin = min(R);
elseif strcmpi(autoOrigin,'off')
    Origin = 0;
else
    error(' ''autoOrigin'' must be ''on'' or ''of'' ')
end

if Origin==0 && strcmpi(Rscale,'log')
    warning(' The origin cannot be set to 0 if R is expressed on a logarithmic axis. The value ''Rmin'' is used instead')
    Origin = min(R);
end

if isempty(circPos)
    if ~isempty(RtickLabel)
        if numel(RtickLabel)~=Ncircles 
            error(' The radial ticklabel must be equal to Ncircles');
        end
        if any(cellfun(@ischar,RtickLabel)==0)
            error(' The radial ticklabel must be a cell array of characters');
        end
    end   
end

if ~isempty(circPos)
    circPos = unique([min(R),circPos,max(R)]);
end
%% Preliminary checks
% case where dimension is reversed
Nrr = numel(R);
Noo = numel(theta);
if isequal(size(Z),[Noo,Nrr]) && Noo~=Nrr,
    Z=Z';
end

% case where dimension of Z is not compatible with theta and R
if ~isequal(size(Z),[Nrr,Noo])
    fprintf('\n')
    fprintf([ 'Size of Z is : [',num2str(size(Z)),'] \n']);
    fprintf([ 'Size of R is : [',num2str(size(R)),'] \n']);
    fprintf([ 'Size of theta is : [',num2str(size(theta)),'] \n\n']);
    error(' dimension of Z does not agree with dimension of R and Theta')
end
%% data plot
rMin = min(R);
rMax = max(R);
thetaMin=min(theta);
thetaMax =max(theta);
    
% Definition of the mesh
cax = newplot;
Rrange = rMax - rMin; % get the range for the radius
[rNorm] = getRnorm(Rscale,Origin,R,Rrange); % getRnorm is a nested function
YY = (rNorm)'*cosd(theta);
XX = (rNorm)'*sind(theta);

if colBar==1,
    h = contourf(XX,YY,Z,Contours,'parent',cax,'LineColor','none');
elseif colBar==2 || colBar==3,
    h = contourf(XX,YY,log10(Z),log10(Contours),'parent',cax,'LineColor','none');
end

if ~isempty(ncolor)
    cmap = feval(myColorMap,ncolor);
    colormap(gca,cmap);
else
    colormap(gca,myColorMap);
end

% disp([max(R/Rrange),max(rNorm)])

set(cax,'dataaspectratio',[1 1 1]); set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize', fontsize); axis off;
if ~ishold(cax);
    % make a radial grid
    hold(cax,'on')
    % Draw circles and spokes
    createSpokes(thetaMin,thetaMax,Ncircles,circPos,Nspokes);
    createCircles(rMin,rMax,thetaMin,thetaMax,Ncircles,circPos,Nspokes)
end

%% PLot colorbar if specified
if colBar==1,
    c =colorbar('location','WestOutside', 'fontsize', fontsize);
    caxis([quantile(Z(:),0.01),quantile(Z(:),0.99)])
elseif colBar==2, 
    c = colorbar('YTick',log10(Contours),'YTickLabel',Contours,'location', 'WestOutside', 'FontSize', fontsize);
    caxis(log10([Contours(1) Contours(length(Contours))]));
elseif colBar==3, 
    c = colorbar('YTick',log10(Contours),'YTickLabel',Contours,'location', 'WestOutside', 'FontSize', fontsize);
    caxis(log10([Contours(2) Contours(length(Contours))]));
else
    c = [];
end

%% Outputs
nargoutchk(0,2)
if nargout==1,
    varargout{1}=h;
elseif nargout==2,
    varargout{1}=h;
    varargout{2}=c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function createSpokes(thetaMin,thetaMax,Ncircles,circlesPos,Nspokes)
        
        spokeMesh = linspace(thetaMin,thetaMax,Nspokes);
        if isempty(circlesPos)
            circleMesh = linspace(rMin,rMax,Ncircles);
        else
            circleMesh  =  circlesPos;
        end
        contourD = abs((circleMesh - circleMesh(1))/Rrange+R(1)/Rrange);
        
        if strcmpi(typeRose,'meteo')
            cost = cosd(90-spokeMesh); % the zero angle is aligned with North
            sint = sind(90-spokeMesh); % the zero angle is aligned with North
        elseif strcmpi(typeRose,'default')
            cost = cosd(spokeMesh); % the zero angle is aligned with east
            sint = sind(spokeMesh); % the zero angle is aligned with east
        else
            error('"type" must be "meteo" or "default" ');
        end
        
        for kk = 1:Nspokes
            
            X = cost(kk)*(contourD+0.01);
            Y = sint(kk)*(contourD+0.01);
            
            if  Origin==0
                X(1)=Origin;
                Y(1)=Origin;
            end
            plot(X,Y,LineStyle,'color',LineColor,'linewidth',LineWidth,...
                'handlevisibility','off');
            % plot graduations of angles
            % avoid superimposition of 0 and 360
            if and(thetaMin==0,thetaMax == 360),
                if spokeMesh(kk)<360,
                    
                    text(1.05.*contourD(end).*cost(kk),...
                        1.05.*contourD(end).*sint(kk),...
                        [num2str(spokeMesh(kk),3)],...
                        'fontsize', fontsize,'Interpreter','Latex',...
                        'horiz', 'center', 'vert', 'middle');
                end
            else
                text(1.05.*contourD(end).*cost(kk),...
                    1.05.*contourD(end).*sint(kk),...
                    [num2str(spokeMesh(kk),3), '$^\circ$'],...
                    'fontsize', fontsize,'Interpreter','Latex',...
                    'horiz', 'center', 'vert', 'middle'); 
            end
            
        end
    end
    function createCircles(rMin,rMax,thetaMin,thetaMax,Ncircles,circlePos,Nspokes)
        
        if isempty(circlePos)
            if Origin ==0 % if the origin is set at rMin
                contourD = linspace(0,1+R(1)/Rrange,Ncircles);
            else % if the origin is automatically centered at 0
                contourD = linspace(0,1,Ncircles)+R(1)/Rrange;
            end
        else
            
            contourD = circlePos-circlePos(1);
            contourD = contourD./max(contourD)*max(R/Rrange);
            contourD =[contourD(1:end-1)./contourD(end),1]+R(1)/Rrange;
        end
        
        if isempty(circlePos)
            if strcmpi(Rscale,'linear')||strcmpi(Rscale,'lin'),
                tickMesh = linspace(rMin,rMax,Ncircles);
            elseif strcmpi(Rscale,'log')||strcmpi(Rscale,'logarithmic'),
                tickMesh  = logspace(log10(rMin),log10(rMax),Ncircles);
            else
                error('''Rscale'' must be ''log'' or ''linear'' ');
            end
        else
            tickMesh  = circlePos;
            Ncircles = numel(tickMesh);
        end
        

        % define the grid in polar coordinates
        
        
        if strcmpi(typeRose,'meteo')
            angleGrid = linspace(90-thetaMin,90-thetaMax,100);
        elseif strcmpi(typeRose,'default')
            angleGrid = linspace(thetaMin,thetaMax,100);
        else
            error('"type" must be "meteo" or "default" ');
        end
        
        
        
        
        xGrid = cosd(angleGrid);
        yGrid = sind(angleGrid);
        
        spokeMesh = linspace(thetaMin,thetaMax,Nspokes);
        
        % plot circles
        for kk=1:length(contourD)
            X = xGrid*contourD(kk);
            Y = yGrid*contourD(kk);
            plot(X,Y,LineStyle,'color',LineColor,'linewidth',LineWidth);
        end
        % radius tick label
        for kk=1:Ncircles
            
            %Original: position = 0.51.*(spokeMesh(min(Nspokes,round(Ncircles/2)))+...
            %              spokeMesh(min(Nspokes,1+round(Ncircles/2))));
            position = 0.*(spokeMesh(min(Nspokes,round(Ncircles/2)))+...
                spokeMesh(min(Nspokes,1+round(Ncircles/2))));
        
            if strcmpi(typeRose,'meteo'),position = 90-position; end
            
            
            
            if isempty(RtickLabel),
                rtick = {''}; num2str(tickMesh(kk),2);
            else
                rtick = RtickLabel(kk);
            end
            if abs(round(position)) ==90,
                % radial graduations
                text((contourD(kk)).*cosd(position),...
                    (-0.075+contourD(kk)).*sind(position-4),...
                    rtick,'verticalalignment','BaseLine',...
                    'horizontalAlignment', 'center',...
                    'color', TextColor, ...
                    'fontsize', fontsize,'Interpreter','Latex',...
                    'handlevisibility','off','parent',cax); 
                
                % annotate spokes
                text(contourD(end).*0.6.*cosd(position),...
                    0.07+contourD(end).*0.6.*sind(position),...
                    [labelR],'verticalalignment','bottom',...
                    'horizontalAlignment', 'right',...
                    'fontsize', fontsize,'Interpreter','Latex',...
                    'handlevisibility','off','parent',cax,...
                    'color', TextColor);
            else
                % radial graduations
               text((contourD(kk)).*cosd(position),...
                    (contourD(kk)).*sind(position),...
                    rtick,'verticalalignment','BaseLine',...
                    'horizontalAlignment', 'right',...
                    'color', TextColor, ...
                    'fontsize', fontsize,'Interpreter','Latex',...
                    'handlevisibility','off','parent',cax);

                % annotate spokes
                text(contourD(end).*0.6.*cosd(position),...
                    contourD(end).*0.6.*sind(position),...
                    [labelR],'verticalalignment','bottom',...
                    'horizontalAlignment', 'right',...
                    'fontsize', fontsize,'Interpreter','Latex',...
                    'handlevisibility','off','parent',cax,...
                    'color', TextColor);
            end
        end
        
    end
    function [rNorm] = getRnorm(Rscale,Origin,R,Rrange)
        if strcmpi(Rscale,'linear')||strcmpi(Rscale,'lin')
                rNorm = R-R(1)+Origin;
                rNorm = (rNorm)/max(rNorm)*max(R/Rrange);
        elseif strcmpi(Rscale,'log')||strcmpi(Rscale,'logarithmic')
            if rMin<=0
                error(' The radial vector cannot be lower or equal to 0 if the logarithmic scale is used');
            end
            rNorm = log10(R); %normalized radius [0,1]
            rNorm =rNorm-rNorm(1);
            rNorm = (rNorm)/max(rNorm)*max(R/Rrange);
        else
            error('''Rscale'' must be ''log'' or ''linear'' ');
        end
    end
end

