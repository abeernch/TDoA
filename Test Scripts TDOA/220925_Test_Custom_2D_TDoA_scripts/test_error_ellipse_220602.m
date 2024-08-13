clear;clc;close all
C = [1 3.5; 0 3];
[v,lam] = eig(C);
eig = lam(lam~=0);
eig_max = max(eig);
eig_min = min(eig);
[sorted_v,i] = sort(v,'descend');
v_max = v(:,i(2));

gamma = 4.601;
maj_axis = sqrt(gamma*eig_max);
min_axis = sqrt(gamma*eig_min);

theta = linspace(0,2*pi,100);

rot_ang = atan2(v_max(2),v_max(1));

pos_x = maj_axis.*cos(theta);
pos_y = min_axis.*sin(theta);

x = cos(rot_ang)*pos_x - sin(rot_ang)*pos_y;
y = sin(rot_ang)*pos_x + cos(rot_ang)*pos_x;

ellipse = [x;y];
figure; plot(ellipse(1,:),ellipse(2,:))