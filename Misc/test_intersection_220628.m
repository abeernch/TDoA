% x=-pi:pi/100:pi;
% y1=sin(10*x);
% y2=sin(x);
x1 = xiso{2}(1,:);
x2 = xiso{3}(1,:);
y1 = xiso{2}(2,:);
y2 = xiso{3}(2,:);

dy1y2 = y2-y1;                                                  % Subtract To Create Zero-Crossings
dx1x2 = x2-x1;                                                  % Subtract To Create Zero-Crossings
idx_Y = find(diff(sign(dy1y2)));                                  % Find Approximate Indices Of Zero-Crossings
idx_X = find(diff(sign(dx1x2)));                                  % Find Approximate Indices Of Zero-Crossings

for k = 1:numel(idx_X)
    idxrng_X = max(idx_X(k)-2,1) : min(numel(x1),idx_X(k)+2);          % Restrict Index Range To Meet Requirements For _interp1X
    xv(k) = interp1(dy1y2(idxrng_X), idxrng_X, 0);               % Calculate ‘x’ Coordinates Of Zero-Crossings
    
end

for k = 1:numel(idx_Y)
    idxrng_Y = max(idx_Y(k)-2,1) : min(numel(y1),idx_Y(k)+2);          % Restrict Index Range To Meet Requirements For _interp1X
    yv(k) = interp1(idxrng_Y, idxrng_Y, 0);              % Calculate ‘x’ Coordinates Of Zero-Crossings (Should Be Close To Zero For Most, Although Not All)
end

% figure(1)
% plot(x,y1,'b');
% hold on
% plot(x,y2,'g');
% hold off
% grid
% ylim([-2 2])
% 
% figure(2)
% plot(x,y1,'b');
% hold on
% plot(x,y2,'g');
% plot(xv, yv, 'xr')
% hold off
% grid
% axis([-10 1    -1.5  1.5])
% text(xv(end-2:end), yv(end-2:end), compose('\\leftarrow (%7.4f, %7.4f)', [xv(end-2:end); yv(end-2:end)].'), 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Rotation',-15)
% xlabel('x')
% ylabel('y')
% title('Last 3 Intersections With Coordinate Values')
% legend('y_1','y_2','Intersections','Location','NW')