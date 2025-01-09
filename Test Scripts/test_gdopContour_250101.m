function gdop = gdopContour(xmin, xmax, ymin, ymax, SVs)
%gdop = gdopContour(xmin, xmax, ymin, ymax, SVs)
%gdopContour shows a psuedocolor image of the PDOP (GDOP in 2-D is referred
%to as PDOP). xmin/xmax/ymin/ymax are the limits on the plot. SVs is an
%[x,y] matrix showing the position of the "Space Vehicles" (GPS
%terminology). You can have unlimited SVs. 
%Reference: P. Dana at
%http://www.colorado.edu/geography/gcraft/notes/gps/gif/gdop.gif
if(xmax > xmin)
    if(ymax > ymin)
        gdop = [];
        for i = xmin:1000:xmax
            val = [];
            for j = ymin:1000:ymax
                val = vertcat(val, gdop1(i,j, SVs));
            end
            gdop = horzcat(gdop, val);
        end
    end
end
x = [xmin:1000:xmax];
y = [ymin:1000:ymax];
figure, h = pcolor(x,y,gdop);
% shading interp
% set(h, 'CDataMapping', 'Direct')
% colorbar
hold on, scatter(SVs(:,1), SVs(:,2), 'd', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
contour(x,y,gdop)
end

function GDOP = gdop1(x,y, SVs)
%gdop1 is a measure of precision. References:
%http://www.colorado.edu/geography/gcraft/notes/gps/gps.html
%http://www.colorado.edu/geography/gcraft/notes/gps/gif/gdop.gif
pos = [x,y];    %a given x-y coord
%ranges for receiver position estimate
Ri = sqrt((SVs(:,1)-pos(1,1)).^2 + (SVs(:,2)-pos(1,2)).^2);
Dx = (SVs(:,1)-pos(1))./Ri; % directional derivative
Dy = (SVs(:,2)-pos(2))./Ri; % directional derivative
A = horzcat(Dx, Dy); %[Dx0, Dy0; Dx1, Dy1]
P = inv(A'*A);
GDOP = sqrt(trace(P));
end