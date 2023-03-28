function [alphaW] = alphaWLowe2005(wave, vege, canopy)
% Reduction factor alphaW for within canopy velocity u = alphaW * u0
% (equation 2 in Zhu et al., 2022), where u0 is velocity outside the
% canopy. The calculation formula is referred to Lowe et al. (2005).

%% input variables
% wave.h:   [m] water depth;
% wave.Tw:  [s] wave period;
% wave.H0:  [m] incident wave height;
%
% vege.b:   [m] blade width projected to the direction of wave propagation;
% vege.Ac:  [m^2] cross section area of blade;
% vege.d2:  [m] blade length (vegetation height).
% vege.d3:  [m] gap below the canopy with 0 for submerged vegetation and >0
%              for suspended kelp.
% vege.tip: 0 --> fixed at bottom; 1 --> fixed at top
%
% canopy.N: [blades/m^2] canopy density, number of blades per square meter;
%
% functions used:
% [Cd] = CdLN2016(Um, Tw, b)
% [Cm] = CmLN2016(Um, Tw, b)

%% output variables
% alphaW: [-] alphaW in Zhu et al. (2022), the ratio of the within canopy
%            velocity to the above-canopy velocity

%%
h=wave.h;
Tw=wave.Tw;
H = wave.H0;
b=vege.b;
Ac=vege.Ac;
N=canopy.N;
l=vege.d2;
if vege.tip > 0
    z = vege.d2 + vege.d3;
else
    z = vege.d3;
end

k = waveNum(h,Tw);
omega = 2*pi/Tw;
Umax = H/2*omega*cosh(k*z)/sinh(k*h);  %% kz << 1 so cosh(kz)=1;

Cd = CdLN2016(Umax, Tw, b);
Cm = CmLN2016(Umax, Tw, b);
Cf = 0.01;

lambda_f = l*b *N;
lambda_p = Ac *N;

Ld = 2*l*(1-lambda_p)/(Cd*lambda_f);
Ls = 2*l/Cf;
    
U_infw = @(t) Umax*cos(omega * t);
dU_infw_dt = @(t) -omega*Umax*sin(omega * t);

dt = 0.01;
uim1 = 0;
ui= 0;
t = (2*dt:dt:20*Tw)';
nt = length(t);

uw = zeros(nt,1);
uw_inf = zeros(nt,1);
for it = 1:nt
    uip1 = (dU_infw_dt(t(it)) + abs(U_infw(t(it)))*(U_infw(t(it)))/Ls...
        - abs(ui)*(ui)/Ld )....
        * dt/(1+Cm*lambda_p/(1-lambda_p)) -uim1+2*ui;
    ui = uip1;
    uim1 = ui;
    uw(it) = ui;
    uw_inf(it) = U_infw(t(it));
end

alphaW = rms(uw(nt- floor(Tw/dt):nt))/rms(uw_inf(nt- floor(Tw/dt):nt));
