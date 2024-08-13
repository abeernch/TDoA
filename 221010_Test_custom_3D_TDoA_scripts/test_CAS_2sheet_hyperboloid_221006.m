%% Implicit evaluation of 2-sheet hyperboloid and CAS verification of derivation 
clc;clear;
% Shaping parameters of the quadric surface
a = 1; 
b = 1; 
c = 1;

f = @(x,y,z) + (x.^2)./a^2 + (y.^2)./b^2 - (z.^2)./c^2 + 1;

% plotting interval
interval = [-10 10 -10 10 -10 10];

% plot figure
figure
a = fimplicit3(f,interval);
shading interp
%%
clc;clear
syms x y z a b c
eqn = (x.^2)/a^2 + (y^2)/b^2 - (z^2)/c^2 == -1;

s = solve(eqn,z);
s_me = [(b/(a*c))*sqrt(a^2*z.^2 + c^2*x.^2 + a^2*c^2);
        -(b/(a*c))*sqrt(a^2*z.^2 + c^2*x.^2 + a^2*c^2);];

chk = s-s_me;
simplify(expand(chk))

pretty(s)