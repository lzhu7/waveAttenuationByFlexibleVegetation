function [vege, canopy, hydroCoeff]....
    = regularWaveDecayRigidVegetation(wave, vege, canopy, hydroCoeff)

%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2022%%%%%%%%%%%%%%%%%%%%%%
%---------------------------version 2022.04--------------------------------
%
% regularWaveDecayRigidVegetation is a MATLAB transcript to calculate the
% wave attenuation by rigid vegetation (and suspended kelp) in regular
% monochromatic waves. This code is based on the analytical model in Zhu
% and Zou (2017), which was validated by the laboratory experiments in
% Ozeren et al. (2014) for submerged vegetation. Note that the wave decay
% coefficient kD in this code is kD = kv/H0, where kv is the exponential
% decay coefficient defined in Zhu and Zou (2017) and H0 is the incident
% wave height. The wave decay coefficient is applicable for both
% exponential wave decay (Kobayashi et al., 1993) and fractional wave decay
% (Dalrymple et al, 1984) forms. These two wave decay forms are linked
% through a piecewise method in Appendix C of Zhu et al. (2020).
%
% Recommended Citation:
%
% [1] Zhu, L., & Zou, Q. (2017). Three-layer analytical solution for wave
% attenuation by suspended and nonsuspended vegetation canopy. Coast. Eng.
% Proc. 1 (35), 27. https://doi.org/ 10.9753/icce.v35.waves.27.
%
% [2] Zhu, L., Huguenard, K., Fredriksson, D. W., & Lei, J. (2022). Wave
% attenuation by flexible vegetation (and suspended kelp) with blade
% motion: Analytical solutions. Advances in Water Resources, 162, 104148.
% https://doi.org/10.1016/j.advwatres.2022.104148
%
% The analytical solutions for flexible (and rigid) vegetation (and
% suspended kelp) to attenuate regular monochromatic waves and
% frequency-dependent random waves are referred to Zhu et al. (2022) and
% Zhu et al. (2020a), respectively. More precise numerical models for the
% wave attenuation by flexible vegetation/kelp are referred to Zhu et al.
% (2021) with resolving the (high-order) asymmetric motion of
% vegetation/kelp (Zhu et al., 2020b). These models and theories are
% summarized in the PhD dissertation of Zhu (2020).
%
% References:
%
% [1] Dalrymple, R.A., Kirby, J.T., & Hwang, P.A. (1984). Wave diffraction
% due to areas of energy dissipation. J. Waterw. Port Coast. Ocean Eng. 110
% (1), 67–79. https://doi.org/ 10.1061/(ASCE)0733-950X(1984)110:1(67).
%
% [2] Kobayashi, N., Raichle, A.W., & Asano, T. (1993). Wave attenuation by
% vegetation. J. Waterw. Port Coast. Ocean Eng. 119 (1), 30–48.
% https://doi.org/10.1061/(ASCE) 0733-950X(1993)119:1(30).
%
% [3] Ozeren, Y., Wren, D.G., & Wu, W. (2014). Experimental investigation
% of wave attenuation through model and live vegetation. Journal of
% Waterway, Port, Coastal, and Ocean Engineering, 140(5), p.04014019.
%
% [4] Zhu, L. (2020). Wave Attenuation Capacity of Suspended Aquaculture
% Structures with Sugar Kelp and Mussels. PhD Thesis. The University of
% Maine. https://digitalcommons.library.umaine.edu/etd/3222
%
% [5] Zhu, L., Huguenard, K., Zou, Q., Fredriksson, D.W., & Xie, D.
% (2020a). Aquaculture farms as nature-based coastal protection: random
% wave attenuation by suspended and submerged canopies. Coast. Eng. 160,
% 103737 https://doi.org/10.1016/j.coastaleng.2020.103737
%
% [6] Zhu, L., Zou, Q., Huguenard, K., & Fredriksson, D. W. (2020b).
% Mechanisms for the Asymmetric Motion of Submerged Aquatic Vegetation in
% Waves: A Consistent-Mass Cable Model. JGR: Oceans, 125(2).
% https://doi.org/10.1029/2019JC015517
%
% [7] Zhu, L., Lei, J., Huguenard, K., & Fredriksson, D.W. (2021). Wave
% attenuation by suspended canopies with cultivated kelp (Saccharina
% latissima). Coast. Eng., 103947
% https://doi.org/10.1016/j.coastaleng.2021.103947
%
% [8] Zhu, L., Huguenard, K., Fredriksson, D. W., & Lei, J. (2022). Wave
% attenuation by flexible vegetation (and suspended kelp) with blade
% motion: Analytical solutions. Advances in Water Resources, 162, 104148.
% https://doi.org/10.1016/j.advwatres.2022.104148
%
% Contact longhuan.zhu@maine.edu or lzhu7@mtu.edu in case of any questions.
%
%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2022%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%% input variables in metric units (SI) and example
% % wave parameters
%     wave.h = 0.4;      % [m] water depth
%     wave.rho = 1000;   % [k/m^3] water mass density
%     wave.Tw = 2;       % [s] wave period
%     wave.H0 = 0.2;     % [m] wave height
% 
% % blade parameters
%     vege.b = 0.01;     % [m] projected width of vegetation cross-section
%                        % in the direction of wave propagation
%     vege.Ac = 1e-5;    % [m^2] area of vegetation cross section
%     vege.d2 = 0.15;    % [m] vegetation length (or height)
%     vege.d3 = 0;       % [m] gap below the canopy with 0 for submerged
%                        % vegetation and >0 for suspended vegetation/kelp
%     vege.tip = 1;      % 0 --> vegetation is fixed at bottom;
%                        % 1 --> vegetation is fixed at top
% 
% % canopy parameters
%     canopy.Lv = 4;       % [m] canopy length in the direction of wave
%                          % propagation
%     canopy.N = 500;      % [baldes/m^2] canopy/plant density, number of
%                          % blades per unit area
%     canopy.alpha_w = -1; % default value -1: calculate alpha_w using the
%                          % formula in Lowe et al. (2005)
%                          % set a positive value if it is known
%                          % e.g., 1 --> no velocity reduction in canopy
%                       
% % hydrodynamic coefficients
%     hydroCoeff.Cd = -1;    % default value -1: calcualte Cd using the
%                            % formula in Luhar and Nepf (2016).
%                            % set a positive value if it is known
%                            % e.g., 1.95 --> drag coefficient is 1.95
%     hydroCoeff.Cm = -1;    % default value -1: calcualte Cm using the
%                            % formula in Luhar and Nepf (2016).
%                            % set a positive value if it is known
%                            % e.g., 1 --> added mass coefficient is 1
%
% % functions used:
% % [k]=waveNum(h,T)
% % [alphaW] = alphaWLowe2005(wave, vege, canopy)
% % [Cd] = CdLN2016(Um, Tw, b)
%
%%%%%%%%%%%% output variables
% vege.Fxrms:        [N] root mean square of the total horizontal force on
%                     vegetation
% hydroCoeff.CDm:     [-] mean drag coefficient 
% 
% canopy.x:      [m] distance along the canopy with 20 intervals;
% canopy.H(nx):       [m] Wave height along the canopy as a function of x
%
% canopy.kD:     [m^-2] wave decay coefficient defined in Zhu et al.
%                            (2022)
% canopy.CD:     [-] Bulk drag coefficient for wave attenuation

%% initialization
% wave parameters
h=wave.h;
Tw = wave.Tw;
omegaW=2*pi/Tw;
H0 = wave.H0;
rho_w=wave.rho;

% vegetation parameters
b=vege.b;
d2=vege.d2;
d3=vege.d3;
d1=h-d2-d3;
dz=d2/2; %two segments, three points
z=(-d1-d2:dz:-d1)';

% canopy parameters
N=canopy.N;
Lv = canopy.Lv;

%%
k=waveNum(h,Tw);
GammaZ=cosh(k*(h+z))/sinh(k*h);

if canopy.alpha_w < 0 
    canopy.alpha_w = alphaWLowe2005(wave, vege, canopy);
end
Hv = H0 * canopy.alpha_w; % Hv used to consider incanopy velocity

Um=Hv/2*omegaW*GammaZ(2); % velocity at middle along blade length
if hydroCoeff.Cd>0
    CD = hydroCoeff.Cd;
else
    CD = CdLN2016(Um, Tw, b);
end

kD = CD*b*N*k/(9*pi)...
    *(9*sinh(k*(d2+d3)) - 9*sinh(k*d3) + sinh(3*k*(d2+d3)) - sinh(3*k*d3))...
    /(sinh(k*h)*(2*k*h+sinh(2*k*h)));

%% output
% used hydrodynamic coefficient
hydroCoeff.Cdm = mean(CD);

% total horizontal force on blade
theta=0:2*pi/100:2*pi;
nTheta = length(theta);
FD = zeros(nTheta,1);
for iTheta = 1:nTheta
FD(iTheta) = trapz(z, 1/2*CD*rho_w*b*(Hv/2*omegaW*GammaZ).^2.....
    .*abs(cos(theta(iTheta)))....
    .*(cos(theta(iTheta))));
end
vege.Fxrms = rms(FD);

% wave decay
canopy.kD = kD;
canopy.x = 0:Lv/20:Lv;
canopy.H = H0./(1 + canopy.kD*H0*canopy.x);

% bulk drag coefficient CD
canopy.CD = CD;

%% plot results
% figure;clf; hold on;
% plot(canopy.x/Lv, canopy.H/H0,'r')   
% axis([0 1 0 1])
% xlabel('Distance x/L_v')
% ylabel('Wave height H/H_0')    
% title('Wave height decay over rigid vegetation canopy')