function [mu,phi]=Cantilever_f_phi(n,s)
% find the first n \mu l of cantilever beam, natural frequency
% n: number of mode
% s: distance along the length
% output
% mu
% phi

fun = @(mul_) 1+cos(mul_).*cosh(mul_);
% % omega=zeros(n,1);
mu=zeros(n,1);
ns=length(s);
phi=zeros(ns,n);
for i=1:n
    mu(i) = fzero(fun,[pi*(i-1), pi*i])/s(end);
% %     omega(i) = mu^2*sqrt(EI/m);        % natural frequency
    phi(:,i)=(cos(mu(i)*s(end))+cosh(mu(i)*s(end)))*(sin(mu(i)*s)-sinh(mu(i)*s))...
        +(sin(mu(i)*s(end))+sinh(mu(i)*s(end)))*(cosh(mu(i)*s)-cos(mu(i)*s));
end
