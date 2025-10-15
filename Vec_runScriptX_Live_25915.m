% Load data file as a struct for testing the time performance
S = load('D:\SC2\Projects\TDOA\MS17 TDoA Integrated Scripts\Matlab_testing\data.mat');
fn = fieldnames(S);
raw = S.(fn{1});

step_s = 60000;
numSamples = min(step_s, size(raw,2));

TOA1 = raw(1,1:numSamples);
TOA2 = raw(2,1:numSamples);
TOA3 = raw(3,1:numSamples);
TOA4 = raw(4,1:numSamples);

tic;

toa_arr = [TOA1;TOA2;TOA3;TOA4];

%% All Baslelines
% Get TDoAs
    tdoaPairs(:,2) = (2:(nSensors))';tdoaPairs(:,1) = 1;
    for n = 1:size(tdoaPairs,1)
        tdoa{n} = (toa_arr(1,:) - toa_arr(n+1,:))./200e6;
        [tdoaClean{n},idx{n},~,~] = hampel(tdoa{n});
    end


if sum(sensorMask) == 4
    sol_x = [ solHyp_x1(c,tdoaClean{1},tdoaClean{2},tdoaClean{3});
              solHyp_x2(c,tdoaClean{1},tdoaClean{2},tdoaClean{3}) ];
    sol_y = [ solHyp_y1(c,tdoaClean{1},tdoaClean{2},tdoaClean{3});
              solHyp_y2(c,tdoaClean{1},tdoaClean{2},tdoaClean{3}) ];
    sol_z = [ solHyp_z1(c,tdoaClean{1},tdoaClean{2},tdoaClean{3});
              solHyp_z2(c,tdoaClean{1},tdoaClean{2},tdoaClean{3}) ];

    [solx_1,soly_1] = local2latlon(real(sol_x(1,:)),real(sol_y(1,:)),0,[lat_node{1},lon_node{1},0]);
    [solx_2,soly_2] = local2latlon(real(sol_x(2,:)),real(sol_y(2,:)),0,[lat_node{1},lon_node{1},0]);

    Geosol_x = [solx_1;solx_2];
    Geosol_y = [soly_1;soly_2];
    Geosol_z = sol_z;

elseif sum(sensorMask) == 3
    sol_x = [ solHyp_x1(c,tdoaClean{1},tdoaClean{2});
              solHyp_x2(c,tdoaClean{1},tdoaClean{2}) ];
    sol_y = [ solHyp_y1(c,tdoaClean{1},tdoaClean{2});
              solHyp_y2(c,tdoaClean{1},tdoaClean{2}) ];
    sol_z = -ones(2,numSamples);

    [solx_1,soly_1] = local2latlon(real(sol_x(1,:)),real(sol_y(1,:)),0,[lat_node{1},lon_node{1},0]);
    [solx_2,soly_2] = local2latlon(real(sol_x(2,:)),real(sol_y(2,:)),0,[lat_node{1},lon_node{1},0]);

    Geosol_x = [solx_1;solx_2];
    Geosol_y = [soly_1;soly_2];
    Geosol_z = sol_z;
else
    Geosol_x=[]; Geosol_y=[]; Geosol_z=[];
end

elapsedTime = toc;
fprintf('Execution time for %d samples: %.3f ms\n', numSamples, elapsedTime*1000);

a = Geosol_x;
b = Geosol_y;
cc = Geosol_z;

save('a.mat','a')
save('b.mat','b')
save('c.mat','cc')
clearvars sol_x sol_y sol_z Geosol_z cc