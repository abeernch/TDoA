function iso = isochrone_func(s1,s2,rdoa,x_lim,points,theta)

a = abs(rdoa)/2;
d = norm(s1-s2)/2;
b = sqrt(d^2 - a^2);
center = ((s2-s1)/2);

%% Create generic Hyperbola 
xmax = x_lim;
x = linspace(-xmax,xmax,points); % array of x values for plot
y = (a/b).*(sqrt(x.^2 + b^2)); % corresponding y values

%% Apply coord transformation
for i = 1:length(x)
    [xout,yout] = coord_tfm(x(i),y(i),theta,center(1),center(2),rdoa);
    x_trf(i) = xout; 
    y_trf(i) = yout;
end

%% Choose hyperbola leaf
if rdoa < 0
    iso = [x_trf;y_trf]+s1;
elseif rdoa >= 0
    iso = [-x_trf;-y_trf]+s1;
end
end