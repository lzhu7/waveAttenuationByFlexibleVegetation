function [k]=waveNum(h,T)
% Solve for wave number based on linear wave theory:
% (2*pi/T)^2=g*k*tanh(k*h)
% h: water depth [m];
% T: wave period [s];
% k: wave number;

%%
g=9.81;
omega=2*pi./T; % angular frequency
f=@(k) g*k*tanh(k*h)-omega^2; % dispersion relationship function
k_deep=omega^2/g; % initial guss based on deep-water anular freq.
k_shallow=omega/sqrt(g*h);
if abs(f(k_shallow))<= abs(f(k_deep))
    k=fzero(f,k_shallow);
else
    k=fzero(f,k_deep);
end

if k<=0
    error(['wave number is ',num2str(k), '<=0 '])
end