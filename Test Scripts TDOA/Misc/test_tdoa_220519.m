clear 
close all
%% TDOA Implementation
%parameters
range_s = 50e3; %sensor range for each dimension
range_T = 100e3; %target range for each dimension
c = 3e8;    %the speed of light
M = 4;      %number of sensors
E = linspace(50,0.1,1e3); 
in_est_error = 0;  %initial estimate error std in meter
trials = 1;
fail_thr = 10e3;   % threshold for fail of estimator
RMSE = zeros(length(E),1);
for m=1:length(E)
rmse = zeros(trials,1);
err_std = E(m)*1e-9; %tdoa measurement error std in sec
for k=1:trials
%sensor positon vectors
P = zeros(3,M); %includes all sensor positon vectors
% for ii=1:M
%     P(:,ii)=range_s*2*(rand(3,1)-0.5);
    P = [0 5e3 10e3 5e3; 0 5e3 0 -5e3;0 1e2 0.5e2 -1e2]; 
% end
 
% p_T = range_T*2*(rand(3,1)-0.5);   %target positon vector
p_T = [10e3;10e3;0.25e2];
%finding TOAs 
dummy = repmat(p_T,1,M)-P;
toa = zeros(M,1);   %includes all toa information
for ii = 1:M
    toa(ii) = norm(dummy(:,ii))/c;    
end
tdoa = toa-toa(1); tdoa(1)=[];
tdoa = tdoa + err_std*randn(M-1,1);

%%% Taylor Series Expansion Solution
p_T_0 = p_T + in_est_error*randn(3,1);    %initial estimate with some error (penalty term)
d = c*tdoa;
f = zeros(M-1,1);
del_f = zeros(M-1,3);
for ii=2:M
   f(ii-1)=norm(p_T_0-P(:,ii))-norm(p_T_0-P(:,1)); 
   del_f(ii-1,1) = (p_T_0(1)-P(1,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(1)-P(1,1))*norm(p_T_0-P(:,1))^-1;
   del_f(ii-1,2) = (p_T_0(2)-P(2,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(2)-P(2,1))*norm(p_T_0-P(:,1))^-1;
   del_f(ii-1,3) = (p_T_0(3)-P(3,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(3)-P(3,1))*norm(p_T_0-P(:,1))^-1;    
end
x_nonlin = pinv(del_f)*(d-f)+p_T_0;
rmse(k) = norm(p_T-x_nonlin)^2;
end
fails = sum(rmse > fail_thr^2);
RMSE(m) = sqrt(mean(rmse(rmse < fail_thr^2)));
end 
figure
plot(E*1e-10,RMSE);
ylabel('RMSE (m)');
xlabel('\sigma_e (sec)')
%%
% 
% %shows the sensor position, target and estimation 
figure
plot3(P(1,:), P(2,:),P(3,:),'o'); hold on;
plot3(p_T(1), p_T(2),p_T(3),'k*');
xlim([-range_T range_T]);ylim([-range_T range_T]);zlim([-range_T range_T]);
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
plot3(x_nonlin(1), x_nonlin(2),x_nonlin(3),'ms','MarkerSize',7.75); 
legend('Sensor Positions', 'Target Position', 'Target Estimation')
grid on; 
hold off;
xlim([-10e3 10e3]); ylim([-10e3 10e3]); zlim([-1e3 1e3])