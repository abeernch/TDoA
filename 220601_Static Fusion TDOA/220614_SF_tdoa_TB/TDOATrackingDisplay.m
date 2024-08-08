classdef TDOATrackingDisplay < matlab.System
    % This is a class to visualize the scenario, detections and
    % tracks in the TDOA tracking.

    % Some properties to allow tuning
    properties (Nontunable)
        PositionSelector = [1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0];
        VelocitySelector = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1];
        XLimits = [-5000 5000];
        YLimits = [-5000 5000];
        LogAccuracy = false;
        PlotGDOP = true;
        Title = '';
    end

    % Plotters (Read-only)
    properties (SetAccess = protected)
        Figure
        Axes
        TheaterPlot
        TrajectoryPlotter
        EmitterPlotter
        ReceiverPlotter
        TrackPlotter
        FusedDetectionPlotter
        TDOAPlotters
    end

    % For plotting track vs detection accuracy for single-object use-case.
    properties (Access = protected)
        LoggedFusedMeasurementAccuracy
        LoggedTrackAccuracy
        LoggedTime
    end

    methods
        function obj = TDOATrackingDisplay(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods (Access = protected)
        function setupImpl(obj, scenario, rxPairs)
            f = figure('Units','normalized','Position',[0.1 0.1 0.8 0.8],'HandleVisibility','on');
            ax = axes(f);
            tp = theaterPlot('Parent',ax,'XLimits',obj.XLimits,'YLimits',obj.YLimits,'AxesUnits',["km","km","km"]);
            obj.Figure = f;
            obj.Axes = ax;
            obj.TheaterPlot = tp;
            clrs = darkColorOrder;
            colororder(ax, clrs)
            obj.EmitterPlotter = platformPlotter(tp,'MarkerFaceColor',clrs(2,:),'Marker','^','DisplayName','Emitters');
            obj.ReceiverPlotter = platformPlotter(tp,'MarkerFaceColor',clrs(3,:),'Marker','^','DisplayName','Receivers');
            obj.FusedDetectionPlotter = detectionPlotter(tp,'MarkerFaceColor',clrs(4,:),'Marker','o','DisplayName','Static Fused Position');
            obj.TrackPlotter = trackPlotter(tp,'MarkerFaceColor',clrs(1,:),'Marker','s','DisplayName','Tracks','ConnectHistory','on','ColorizeHistory','on');
            hLine = findall(gca,'Tag','tpTrackHistory_Tracks');
            hLine.LineWidth = 3;
            hold (ax,'on');
            obj.TDOAPlotters = gobjects(13,1);
            for i = 1:13
                obj.TDOAPlotters(i) = plot(nan,nan,'Color',clrs(i,:),'LineWidth',1);
            end
            trajP = trajectoryPlotter(tp,'Color','k','LineWidth',1,'LineStyle','-');
            r = record(clone(scenario));
            poses = horzcat(r.Poses);
            pos = cell(size(poses,1),1);
            for i = 1:size(poses,1)
                pos{i} = vertcat(poses(i,:).Position);
            end
            trajP.plotTrajectory(pos);
            title(ax,obj.Title);
            if obj.PlotGDOP
                if isvector(rxPairs)
                    pairIds = zeros(0,2);
                    n = numel(rxPairs);
                    for i = 1:n
                        thisN = n - i;
                        thisPair = [i*ones(thisN,1) ((i+1):n)'];
                        pairIds = [pairIds;thisPair]; %#ok<AGROW> 
                    end
                    rxPairs = rxPairs(pairIds);
                end
                x = linspace(obj.XLimits(1),obj.XLimits(2),500);
                y = linspace(obj.YLimits(1),obj.YLimits(2),500);
                invGDOP = computeGDOP(scenario, rxPairs,x,y);
                imagesc(ax, invGDOP,'XData',x,'YData',y);            
                colorbar
            end
        end

        function stepImpl(obj, scenario, tdoaPairs, tdoaDetections, fusedDetections, tracks)
            platIDs = cellfun(@(x)x.PlatformID,scenario.Platforms);
            isReceiver = ismember(platIDs,tdoaPairs(:));
            sIdx = cellfun(@(x)x.SensorIndex,tdoaDetections);
            uqSIdx = unique(sIdx);
            for i = 1:numel(uqSIdx)
                plotTDOADetection(obj, obj.TDOAPlotters(i), tdoaDetections(sIdx == uqSIdx(i)));
            end
            poses = platformPoses(scenario);
            emitPoses = poses(~isReceiver);
            receivePoses = poses(isReceiver);
            obj.EmitterPlotter.plotPlatform(vertcat(emitPoses.Position));
            obj.ReceiverPlotter.plotPlatform(vertcat(receivePoses.Position));
            allDets = [fusedDetections{:}];
            if ~isempty(allDets)
                triangulatedPos = horzcat(allDets.Measurement);
                triangulatedPosCov = cat(3,allDets.MeasurementNoise);
            else
                triangulatedPos = zeros(3,0);
                triangulatedPosCov = zeros(3,3,0);
            end
            obj.FusedDetectionPlotter.plotDetection(triangulatedPos',triangulatedPosCov);
            [pos, posCov] = getTrackPositions(tracks,obj.PositionSelector);
            vel = getTrackVelocities(tracks,obj.VelocitySelector);
            trackIDs = string([tracks.TrackID]);
            obj.TrackPlotter.plotTrack(pos,posCov,vel,trackIDs);
            drawnow;
            if obj.LogAccuracy
                if isscalar(fusedDetections) && isscalar(tracks)
                    obj.LoggedFusedMeasurementAccuracy(end+1) = sqrt(triangulatedPosCov(1,1) + triangulatedPosCov(2,2));
                    obj.LoggedTrackAccuracy(end+1) = sqrt(posCov(1,1) + posCov(2,2));
                    obj.LoggedTime(end+1) = scenario.SimulationTime;
                end
            end
        end

        function plotTDOADetection(obj, plotter, tdoaDetections)
            x = zeros(0,1);
            y = zeros(0,1);
            for i = 1:numel(tdoaDetections)
                [thisX, thisY] = tdoaHyperbola(tdoaDetections{i},'2D');
                x = [x;nan;thisX(:)]; %#ok<AGROW>
                y = [y;nan;thisY(:)]; %#ok<AGROW>
            end
            plotter.XData = x;
            plotter.YData = y;
        end
    end

    methods
        function zoomOnTrack(display, trackID)
            ax = gca;
            xlim = ax.XLim;
            ylim = ax.YLim;
            txt = findall(ax,'Type','Text','String',string(trackID));
            if ~isempty(txt)
                ax.XLim = txt.Position(1) + [-250 250];
                ax.YLim = txt.Position(2) + [-250 250];
                drawnow;
                snap = takesnap(display.Figure);
                showsnap(snap);
            end
            ax.XLim = xlim;
            ax.YLim = ylim;
        end

        function plotDetectionVsTrackAccuracy(display)
            f = figure('Units','normalized','Position',[0.1 0.1 0.6 0.6]);
            ax = axes(f);
            plot(ax,display.LoggedTime,display.LoggedFusedMeasurementAccuracy,'LineWidth',2);
            hold (ax, 'on');
            plot(ax, display.LoggedTime,display.LoggedTrackAccuracy,'LineWidth',2);
            legend(ax, 'Fused Position','Track');
            title(ax,'Target Position Accuracy');
            xlabel(ax,'Time (s)');
            ylabel(ax,'Position Uncertainty (m)');
            grid(ax,'on');
        end
    end

    methods (Static)
        function plotConstantTDOACurves()
            r1 = [-10e3;0;0];
            r2 = [10e3;0;0];
            measParams(1).OriginPosition = r1;
            measParams(2).OriginPosition = r2;
            tdoaMax = norm(r2 - r1)/physconst('lightspeed');
            tdoas = linspace(-0.9*tdoaMax, 0.9*tdoaMax,10);
            tdoas = tdoas*1e9;
            xTotal = zeros(0,1);
            yTotal = zeros(0,1);
            textPos = zeros(numel(tdoas),2);
            for i = 1:numel(tdoas)
                tdoaDet = objectDetection(0,tdoas(i),'MeasurementParameters',measParams);
                [x,y] = tdoaHyperbola(tdoaDet,'3D');
                xTotal = [xTotal;nan;x];
                yTotal = [yTotal;nan;y];
                textPos(i,1) = x(85) - sign(tdoas(i))*2e3;
                textPos(i,2) = y(85) - 1e3;
            end

            f = figure('Units','normalized');
            f.Position = [0.1 0.1 0.8 0.8];
            clrs = lines(2);
            plot(xTotal/1e3,yTotal/1e3,'LineWidth',2);
            hold on;
            plot(r1(1)/1e3,r1(2)/1e3,'o','Color',clrs(2,:));
            plot(r2(1)/1e3,r2(2)/1e3,'o','Color',clrs(2,:));
            text(r1(1)/1e3-2,r1(2)/1e3,'R_1')
            text(r2(1)/1e3+0.8,r2(2)/1e3,'R_2')
            for i = 1:numel(tdoas)
                text(textPos(i,1)/1e3,textPos(i,2)/1e3,['$$',sprintf('%0.2f',tdoas(i)/1e3),'\mu s $$'],'Interpreter','latex','HorizontalAlignment','center','FontWeight','bold','FontSize',10);
            end
            xlim([-30 30]);
            ylim([-30 30]);
            xlabel('X (km)');
            ylabel('Y (km)');
        end

        function plotGDOPAccuracyMap()
            [scenario, rxPairs] = helperCreateSingleTargetTDOAScenario(3);
            x = linspace(-1e4,1e4,200);
            y = linspace(-1e4,1e4,200);
            [invGDOP, receiverPos] = computeGDOP(scenario, rxPairs, x, y);
            
            f = figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
            ax = axes(f);
            imagesc(invGDOP,'XData',x/1e3,'YData',y/1e3,'Parent',ax);
            ticks = flip([30 50 100 1000]);
            h = colorbar('Ticks',1./ticks,...
                'TickLabels',string(ticks));
            h.Label.String = 'Position Uncertainty (m)';
            hold on;
            clrs = lines(2);
            plot(receiverPos(1,:)/1e3,receiverPos(2,:)/1e3,'^','Color',clrs(2,:),...
                'DisplayName','Receivers','MarkerFaceColor',clrs(2,:));
            legend;
            xlabel('X (km)');
            ylabel('Y (km)');
            ax.YDir = 'normal';
            title(ax,'Geometric Dilution of Precision');
        end
    end
end

function varargout = tdoaHyperbola(tdoaDetection,type)
emissionSpeed = physconst('Lightspeed');
timeScale = 1e9;

% Location of receiver in global coordinates
p1 = tdoaDetection.MeasurementParameters(1);
r1Pos = p1.OriginPosition;

% Location of receiver 2 in global coordinates
p2 = tdoaDetection.MeasurementParameters(2);
r2Pos = p2.OriginPosition;

theta = linspace(-pi/2,pi/2,100);
if strcmpi(type,'2D')
    phi = 0;
else
    phi = linspace(-pi,pi,50);
end

[THETA, PHI] = meshgrid(theta,phi);

% Compute hyperbola quantities
D = norm(r2Pos - r1Pos)/2;
c = tdoaDetection.Measurement*emissionSpeed/timeScale;

% Hyperbola equation in a frame where sensors lie on x axis and are located
% at -D and D on x axis
xC = -c./cos(THETA)./2;
yC = sqrt(4*D^2 - c^2).*tan(THETA).*cos(PHI)./2;
zC = sqrt(4*D^2 - c^2).*tan(THETA).*sin(PHI)./2;

% Translate and rotate to the scenario frame
r0 = (r1Pos + r2Pos)/2;
a = [1;0;0];
b = (r1Pos - r2Pos);
b = b/norm(b);
v = cross(a,b);
s = norm(v);
c = dot(a,b);
V = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
if abs(s) > 0
    R = eye(3) + V + V^2*(1 - c)/s^2;
else
    R = eye(3);
end

x = R(1,1).*xC + R(1,2).*yC + R(1,3).*zC;
y = R(2,1).*xC + R(2,2).*yC + R(2,3).*zC;
z = R(3,1).*xC + R(3,2).*yC + R(3,3).*zC;
x = x + r0(1);
y = y + r0(2);
z = z + r0(3);

if strcmpi(type,'2D')
    x = x(:);
    y = y(:);
    z = z(:);
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
else
    varargout{1} = surf2patch(x,y,z);
end

end

function [invGDOP, receiverPos] = computeGDOP(scenario, rxPairs, x, y)
platIDs = cellfun(@(x)x.PlatformID,scenario.Platforms);
isReceiver = ismember(platIDs,rxPairs(:));
receivers = scenario.Platforms(isReceiver);
receiverPos = zeros(3,numel(receivers));
for i = 1:numel(receivers)
    receiverPos(:,i) = pose(receivers{i},'true').Position;
end
[X, Y] = meshgrid(x,y);
Z = zeros(size(X));
targetPos = [X(:) Y(:) Z(:)]';
emissionSpeed = physconst('Lightspeed');
n = size(rxPairs,1);
invcovxx = zeros(1,numel(X));
invcovyy = zeros(1,numel(X));
invcovxy = zeros(1,numel(X));
invGDOP = zeros(size(X));
tdoaAccuracy = 1e4*(1e-18);
for i = 1:n
    dr1 = pose(scenario.Platforms{platIDs == rxPairs(i,1)}).Position;
    dr2 = pose(scenario.Platforms{platIDs == rxPairs(i,2)}).Position;
    dx1 = targetPos(1,:) - dr1(1);
    dy1 = targetPos(2,:) - dr1(2);
    dz1 = targetPos(3,:) - dr1(3);
    dx2 = targetPos(1,:) - dr2(1);
    dy2 = targetPos(2,:) - dr2(2);
    dz2 = targetPos(3,:) - dr2(3);
    r1 = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
    r2 = sqrt(dx2.^2 + dy2.^2 + dz2.^2);
    tx = (dx1./r1 - dx2./r2)./emissionSpeed;
    ty = (dy1./r1 - dy2./r2)./emissionSpeed;
    invcovxx = invcovxx + tx.^2./tdoaAccuracy;
    invcovyy = invcovyy + ty.^2./tdoaAccuracy;
    invcovxy = invcovxy + tx.*ty./tdoaAccuracy;
end
covxx = invcovyy./(invcovxx.*invcovyy - invcovxy.^2);
covyy = invcovxx./(invcovxx.*invcovyy - invcovxy.^2);
invGDOP(:) = 1./(sqrt(covxx + covyy));

end


function snap = takesnap(f)
children = f.Children;
snap = {copy(children)};
end

function h = showsnap(snap)
h = figure("IntegerHandle","off","Units","normalized","Position",[0.1 0.1 0.8 0.8]);
children = snap{1};
copyobj(children,h);
end

function colorOrder = darkColorOrder
colorOrder = [1.0000    1.0000    0.0667
    0.0745    0.6235    1.0000
    1.0000    0.4118    0.1608
    0.3922    0.8314    0.0745
    0.7176    0.2745    1.0000
    0.0588    1.0000    1.0000
    1.0000    0.0745    0.6510
    1.0000    0.0000    0.0000
    0.4000    0.7804    0.0706
    0.7294    0.8745    1.0000];

colorOrder(11,:) = [1 1 1];
colorOrder(12,:) = [0 0 0];
colorOrder(13,:) = 0.7*[1 1 1];
end
