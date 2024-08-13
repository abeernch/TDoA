clear;close all;clc
%% orth offset of hyperbolae
xmax = 100; ymax = 100;
x = linspace(-xmax,xmax,100); % array of x values for plot
a = 5; b = 3;
y = (a/b).*(sqrt(x.^2 + b^2)); % corresponding y values
% figure(1)
% plot(x,y)
% hold on 
% % plot(x,-y) % Plot other half of hyperbola
% axis([-xmax xmax -ymax ymax]) % specify axis limits
% xlabel('x')
% ylabel('y')

%% Apply coord transformation
for i = 1:length(x)
    [xout,yout] = coord_tfm(x(i),y(i),90,1,1);
    x_trf(i) = xout; 
    y_trf(i) = yout;
end
figure;plot(x_trf,y_trf)
axis([-xmax xmax -ymax ymax]) % specify axis limits

% hold on
% plot(-x_trf,-y_trf)
% hold on
% plot(x,y)
% plot(x,-y)