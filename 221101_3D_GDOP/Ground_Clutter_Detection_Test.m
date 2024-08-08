%% PARAMETERS
clear
c = physconst('lightspeed');fc = 1e9;prf = 8e3; 
tgtpos = [5000 49637; 0 14225; 0 8000]; 
tgtvel = [-10000 120; 0 100; 0 0]; 
tgtrcs = 1;
%% RADAR TRANSCEIVER
rdr = radarTransceiver;rdr.TransmitAntenna.OperatingFrequency = fc;
rdr.ReceiveAntenna.OperatingFrequency = fc;
rdr.Transmitter.Gain = 20;rdr.Transmitter.PeakPower = 1200;

rdr.Waveform.PRF = prf;

sampleRate = 2.4e6;
rdr.Receiver.SampleRate = sampleRate;
rdr.Waveform.SampleRate = sampleRate;
rdr.Waveform.PulseWidth = 10/sampleRate;
bw = [30 30];
% Create a custom element from the pattern of a URA
N = round(0.8859*2./(bw(:).'*pi/180));
N = flip(N);
lambda = freq2wavelen(fc,c);
array = phased.URA(N,lambda/2);
az = -180:.4:180;
el = -90:.4:90;
G = pattern(array,fc,az,el,'Type','efield','normalize',false);
M = 20*log10(abs(G));
P = angle(G);
E = phased.CustomAntennaElement('FrequencyVector',[fc-1 fc+1],...
'AzimuthAngles',az,'ElevationAngles',el,'MagnitudePattern',M,'PhasePattern',P);
rdr.TransmitAntenna.Sensor = E;
rdr.ReceiveAntenna.Sensor = E;
clear G M P az el array lambda N E

% Turn off thermal noise
rdr.Receiver.NoiseMethod = 'Noise power';
rdr.Receiver.NoisePower =db2pow(-110);
rdr.Receiver.Gain = 40;
rdr.Receiver.NoiseFigure = 4;
rdr.NumRepetitions = 10;
%% Scenario & Platforms
rdr.MountingAngles = [0 -2 0];
tgtmotion = phased.Platform('InitialPosition',tgtpos(:,1),'Velocity',tgtvel(:,1));
target = phased.RadarTarget('MeanRCS',tgtrcs*ones(1,1),'OperatingFrequency',fc);  

scene = radarScenario('UpdateRate',prf);
radarplat = platform(scene,'Position',[0 0 50],'Sensors',{rdr});
tgtplat = platform(scene,'Trajectory',...
    kinematicTrajectory('SampleRate',prf,'Position',tgtpos(:,1).',...
    'Velocity',tgtvel(:,1).','AccelerationSource','Property',...
    'AngularVelocitySource','Property'),'signatures',rcsSignature('Pattern',0*ones(2,2)));
rx = phased.ReceiverPreamp('Gain',40,'NoiseFigure',7,'EnableInputPort',true);

%% clutter
[x,y] = meshgrid(linspace(-10000,10000,201));
ht1 = 40*exp(-(x.^2 + y.^2)/30^2);
ht2 = 100*exp(-((x-6200).^2 + y.^2)/1000^2);
ht = ht1 + ht2;
% ht = -1000*ones(201,201);
% ht = zeros(201,201);

h11 = 1000*exp(-((x-3200).^2 + y.^2)/1000^2);
h111 = 1000*exp(-((x-3200).^2 + (y-3200).^2)/1000^2);
h112 = 1000*exp(-((x-3200).^2 + (y+3200).^2)/1000^2);
h12 = 40*exp(-((x+1200).^2 + y.^2)/4000^2);
% h12 = 10*exp(-((x+6200).^2 + y.^2)/4000^2);
h13 = 10*exp(-((x).^2 + (y-6200).^2)/4000^2);
h14 = 10*exp(-((x).^2 + (y+6200).^2)/4000^2);
ht = h11+h111+h112+h12+h13+h14;
figure,subplot(211),mesh(x,y,ht),view([0 90])
%%
% refl = surfaceReflectivityLand(Model = 'ConstantGamma');
refl = surfaceReflectivityLand(Model = 'Billingsley',LandType= 'Desert');

srf = landSurface(scene,RadarReflectivity = refl, ...
    Terrain = ht,Boundary = [-10000,10000;-10000,10000]);
%%
% refl = surfaceReflectivityLand('Model','ConstantGamma','Gamma',-10);
% srf = landSurface(scene,'RadarReflectivity',refl);
clutterGenerator(scene,rdr,'Resolution',100,'RangeLimit',18.75e3,'UseBeam',false)
clut = getClutterGenerator(scene,rdr);
ringClutterRegion(clut,0,10000,100,0)
iqsig = receive(scene);
%%
% figure,
subplot(212),mesh(db(abs(iqsig{1,1})))
view([90 0])
zlim([-140 20])