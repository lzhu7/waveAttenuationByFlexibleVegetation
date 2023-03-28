function [vege, canopy, hydroCoeff]...
    = regularWaveDecayFlexibleVegetation....
    (wave,vege,canopy,hydroCoeff,ctlparameter)

%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2022%%%%%%%%%%%%%%%%%%%%%%
%---------------------------version 2022.04--------------------------------
%
% regularWaveDecayFlexibleVegetation is a MATLAB transcript to calculate
% the wave attenuation by flexible vegetation (and suspended kelp) in
% regular monochromatic waves. This code is based on the analytical model
% in Zhu et al. (2022), which was validated by laboratory experiments (104
% cases) in Luhar et al. (2017) and Lei and Nepf (2019) for submerged
% vegetation and Zhu et al. (2021) for suspended kelp.
%
% Recommended Citation:
% Zhu, L., Huguenard, K., Fredriksson, D. W., & Lei, J. (2022). Wave
% attenuation by flexible vegetation (and suspended kelp) with blade
% motion: Analytical solutions. Advances in Water Resources, 162, 104148.
% https://doi.org/10.1016/j.advwatres.2022.104148
%
% The analytical solutions for frequency-dependent random wave attenuation
% by flexible vegetation can be found in Zhu et al. (2020a). More precise
% numerical models for the wave attenuation by flexible vegetation/kelp can
% be found in Zhu et al. (2021) with resolving the (high-order) asymmetric
% motion of vegetation/kelp (Zhu et al., 2020b). These models and theories
% are summarized in the PhD dissertation of Zhu (2020).
%
% References:
%
% [1] Lei, J., & Nepf, H. (2019). Wave damping by flexible vegetation:
% connecting individual vege dynamics to the meadow scale. Coast. Eng.
% 147, 138–148. https://doi.org/10.1016/j.coastaleng.2019.01.008.
%
% [2] Luhar, M., Infantes, E., & Nepf, H. (2017). Seagrass vege motion
% under waves and its impact on wave decay. J. Geophys. Res.: Oceans 122
% (5), 3736–3752. https://doi.org/10.1002/2017JC012731
%
% [3] Zhu, L. (2020). Wave Attenuation Capacity of Suspended Aquaculture
% Structures with Sugar Kelp and Mussels. PhD Thesis. The University of
% Maine. https://digitalcommons.library.umaine.edu/etd/3222
%
% [4] Zhu, L., Huguenard, K., Zou, Q., Fredriksson, D.W., & Xie, D.
% (2020a). Aquaculture farms as nature-based coastal protection: random
% wave attenuation by suspended and submerged canopies. Coast. Eng. 160,
% 103737 https://doi.org/10.1016/j.coastaleng.2020.103737
%
% [5] Zhu, L., Zou, Q., Huguenard, K., & Fredriksson, D. W. (2020b).
% Mechanisms for the Asymmetric Motion of Submerged Aquatic Vegetation in
% Waves: A Consistent-Mass Cable Model. JGR: Oceans, 125(2).
% https://doi.org/10.1029/2019JC015517
%
% [6] Zhu, L., Lei, J., Huguenard, K., & Fredriksson, D.W. (2021). Wave
% attenuation by suspended canopies with cultivated kelp (Saccharina
% latissima). Coast. Eng., 103947
% https://doi.org/10.1016/j.coastaleng.2021.103947
%
% Contact longhuan.zhu@maine.edu or lzhu7@mtu.edu in case of any questions.
%
%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2022%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%% input variables in metric units (SI) and examples
% % wave parameters
%     wave.h = 0.16;       % [m] water depth
%     wave.rho = 1000;     % [k/m^3] water mass density
%     wave.Tw = 1.4286;    % [s] wave period
%     wave.H0 = 0.028;     % [m] wave height
% 
% % vegetation parameters
%     vege.rho = 920;    % [kg/m^3] mass density of vegetation (individual)
%     vege.EI = 7.5e-8;  % [Nm^2] flexural rigidity of vegetation
%                        % EI = E * I, where E is the flexural modulus [Pa]
%                        % and I is the moment of inertia (second moment of
%                        % area [m^4]. If the material exhibits Isotropic
%                        % behavior then the flexural modulus is equal to
%                        % the elastic modulus (Young's Modulus).
%                        % The second moment of area for rectangular cross-
%                        % section is I = bt^3/12, where b is blade
%                        % width and t is thickness. For circular
%                        % cross-section, I = pi/64(D^4-d^4), where D is
%                        % outer diameter and d is inner diameter. 
%     vege.b = 0.003;    % [m] projected width of vegetation cross-section
%                        % in the direction of wave propagation
%     vege.Ac = 3e-7;    % [m^2] area of vegetation cross section
%     vege.d2 = 0.13;    % [m] vegetation length (or height).
%     vege.d3 = 0.01;    % [m] gap below the canopy with 0 for submerged
%                        % >0 for suspended.
%     vege.tip = 0;      % 0 --> vegetation is fixed at bottom;
%                        % 1 --> vegetation is fixed at top
% 
% % canopy parameters
%     canopy.Lv = 5;         % [m] canopy length in the direction of wave
%                            % propagation
%     canopy.N = 7200;       % [baldes/m^2] canopy/plant density, number of
%                            % blades per unit area
%     canopy.alpha_eps = 1;  % default value 1: without sheltering.
%                            % For dense kelp in Zhu et al. (2021), the
%                            % sheltering factor ranges from 0.360 to 1.047
%                            % with the mean value of 0.630.
%     canopy.alpha_w = -1;   % default value -1: calculate alpha_w the
%                            % formula in Lowe et al. (2005)
%                            % set a positive value if it is known
%                            % e.g., 1 --> no velocity reduction in canopy
%     canopy.alphaM4kD = -1; % default value -1: calculate alpha_M using
%                            % (32) and (33) in Zhu et al. (2022) for
%                            % submerged vegetation and suspended
%                            % vegetation, respectively.
%                            % set a positive value if it is known
%                            % e.g., 1 --> no modification 
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
% % control parameters (default)
%     ctlparameter.tol = 1e-3;       % tolorence
%     ctlparameter.maxNIter = 1000;  % max iteration number
%     ctlparameter.alpha = 0.85;     % relaxation factor 
%     ctlparameter.nMode = 20;       % number of modes 
%     ctlparameter.nSeg = 20;        % segments devided by mode number
% 
% % functions used:
% % [k]=waveNum(h,T)
% % [alphaW] = alphaWLowe2005(wave, vege, canopy)
% % [Cd] = CdLN2016(Um, Tw, b)
% % [Cm] = CmLN2016(Um, Tw, b)
% % [mu,phi]=Cantilever_f_phi(n,s)
%
%%%%%%%%%% output variables
% vege.t(nt):          [s] time series in one wave period with time interval
%                          of 0.01 s.
% vege.x(nz, nt):      [m] time-dependent horizontal deflection of
%                          vegetation along vertical coordiante;
% vege.z(nz, nt):      [m] time-dependent vertical coorindate of vege
%                          segments;
% vege.xAmp(nz):       [m] amplitude of vegetation motion along the length
% vege.Fxrms:          [N] root mean square of the total horizontal force
%                          on vegetation.
% vege.omegaN(nMode):  [Hz] first n natural frequncies of vegetation, which
%                          were presented as lambda_n in Zhu et al. (2022)
%
% hydroCoeff.CDm:      [-] mean drag coefficient; 
% hydroCoeff.CMm:      [-] mean added mass coefficient;
% hydroCoeff.CdFixEnd: [-] drag coefficient at the fixed end of vegetation
% hydroCoeff.CmFixEnd: [-] added mass coefficient at the fixed end of
%                          vegetation
% 
% canopy.x:            [m] distance along the canopy with 20 intervals;
% canopy.H(nx):        [m] Wave height along the canopy as a function of x,
%                          where the nonlinearity of blade motion is not
%                          considered 
% canopy.HModified:    [m] Modified wave height along the canopy by
%                          considering nonlinearity of blade motion
% canopy.kD:           [m^-2] wave decay coefficient defined in Zhu et al. 
%                             (2022)
% canopy.kDModified:   [m^-2] modified wave decay coefficient by considering
%                           the nonlinear effects of blade motion
% canopy.CD:           [-] Bulk drag coefficient for wave attenuation by
%                          vegetation with linear motion
% canopy.CDModified:   [-] Modified bulk drag coefficient for wave
%                          attenuation by vegetation consdiering effects of
%                          nonlinear vegetation motion
% canopy.Le2L:         [-] Effective vege length ratio Le2L = le/l, hwere
%                          le is effective blade length and l is the
%                          original blade length
% canopy.Le2LModified: [-] Modified effective vege length ratio Le2L = le/l
%                          by considering nonlinear effects of blade motion

%% initialization
% control parameters
if nargin < 5
    tol = 1e-3;
    maxNIter = 1000;
    alpha = 0.85; % relaxation factor
    nMode = 20; % >=5
    nSeg = 20; % >=5
else
    tol = ctlparameter.tol;
    maxNIter = ctlparameter.maxNIter;
    alpha = ctlparameter.alpha; %% relaxation factor
    nMode = ctlparameter.nMode;
    nSeg = ctlparameter.nSeg;
end

% wave parameters
h = wave.h;
rho_w = wave.rho;
Tw = wave.Tw;
omegaW = 2*pi/Tw;
H0 = wave.H0;
k = waveNum(h,2*pi/omegaW);

% vegetation parameters
rho_v = vege.rho;
EI = vege.EI;
b = vege.b;
Ac = vege.Ac;
d2 = vege.d2;
d3 = vege.d3;
d1 = h-d2-d3;
dz = d2/(nSeg*nMode);
z = (-d1-d2:dz:-d1)';
nz=length(z);
s = (0:dz:d2)';
omegaN = zeros(nMode,1);
zetaN = zeros(nMode,1);
[muN, phiN]=Cantilever_f_phi(nMode,s);
if vege.tip
    phiN0=phiN;
    for j=1:nMode
        for iz=1:nz
            phiN(iz,j)=phiN0(nz-iz+1,j);
        end
    end
end
GammaZ=cosh(k*(h+z))/sinh(k*h);
gammaS=zeros(nz,1); % transfer functions
gammaC=zeros(nz,1);

% canopy parameters
Lv = canopy.Lv;
N = canopy.N;

%% maine procedure
if canopy.alpha_w < 0 
    canopy.alpha_w = alphaWLowe2005(wave, vege, canopy);
end
Hv = H0 * canopy.alpha_w; % Hv used to consider incanopy velocity
m=(rho_v+rho_w)*Ac; % estimated
mI= m; % estimated
c=rho_w*b*Hv/2*omegaW*8/(3*pi)*trapz(z, GammaZ.^3)/trapz(z, GammaZ.^2);
% estimated
nIter = 0;
Hv0 = Hv/2;
c0 = c/2;
while (abs((Hv0-Hv)/Hv) + norm(c-c0)/norm(c))>tol && nIter<maxNIter
    Hv0 = Hv;
    c0=c;  
    nIter=nIter+1;        
    % transfer function gammaS gammaC
    sumS=zeros(nz,1);
    sumC=zeros(nz,1);
    for j=1:nMode
        omegaN(j)=muN(j)^2 .....
            *sqrt(trapz(s,EI.*phiN(:,j).^2)/trapz(s,m.*phiN(:,j).^2));
        zetaN(j)=(trapz(s,c.*phiN(:,j).^2)/trapz(s,m.*phiN(:,j).^2))...
            ./(2*omegaN(j));
        DN=trapz(z,c.*GammaZ.*phiN(:,j))./trapz(z,m.*phiN(:,j).^2);
        IN=trapz(z,mI.*GammaZ.*phiN(:,j))./trapz(z,m.*phiN(:,j).^2);

        sumS=sumS+phiN(:,j)*(omegaW*IN*(omegaN(j)^2-omegaW^2)...
            -2*zetaN(j)*omegaN(j)*omegaW*DN)...
            /((omegaN(j)^2-omegaW^2)^2....
            +(2*zetaN(j)*omegaN(j)*omegaW)^2);
        sumC=sumC+phiN(:,j)*(DN*(omegaN(j)^2-omegaW^2)...
            +2*zetaN(j)*omegaN(j)*omegaW*omegaW*IN)...
            /((omegaN(j)^2-omegaW^2)^2....
            +(2*zetaN(j)*omegaN(j)*omegaW)^2);
    end
    gammaS = omegaW./GammaZ.*sumS;
    gammaC = omegaW./GammaZ.*sumC;

    % update Cd Cm m mI   
    Um = Hv/2*omegaW*GammaZ.*sqrt((1+gammaS).^2+gammaC.^2);
    if hydroCoeff.Cd>0
        Cd = ones(nz,1).*hydroCoeff.Cd;
    else
        Cd = CdLN2016(Um, Tw, b);
    end
    if hydroCoeff.Cm>0
        Cm = ones(nz,1).*hydroCoeff.Cm;
    else
        Cm = CmLN2016(Um, Tw, b);
    end

    m = rho_v*vege.Ac+Cm*rho_w*Ac;
    mI = (1+Cm)*rho_w*Ac;

    % update linear coefficient c    
    c = 1/2*Cd*rho_w*b*Hv/2*omegaW*8/(3*pi)....
        *trapz(z, GammaZ.^3.*(((1+gammaS).^2+gammaC.^2).^(3/2)))....
        /trapz(z, GammaZ.^2.*((1+gammaS).^2+gammaC.^2))....
        + c0*alpha;
      
    kD = 4*canopy.alpha_eps*b*N*k^2/(3*pi*sinh(k*h)*(2*k*h+sinh(2*k*h)))....
        *trapz(z, Cd.*((1+gammaS).^2+gammaC.^2).^(3/2).*(cosh(k*(h+z))).^3);
    
    Hv = H0/(1 + kD*H0*Lv/2) * canopy.alpha_w; %% Hv at mid canopy
end
disp(['nIter = ',num2str(nIter),': error = ',....
    num2str((abs((Hv0-Hv)/Hv) + norm(c-c0)/norm(c)))]);
if nIter >= maxNIter
    error('DO NOT CONVERGE! Adjust control parameters!')
end

%% output
% output for the used hydrodynamic coefficients
hydroCoeff.Cdm = mean(Cd);
hydroCoeff.Cmm = mean(Cm);
hydroCoeff.CdFixEnd = Cd(end);
hydroCoeff.CmFixEnd = Cm(end);

% vegetation natural frequency
vege.omegaN = omegaN;

% amplitude of vegetation motion
xAmp=Hv/2*GammaZ.*sqrt(gammaS.^2 + gammaC.^2);
vege.xAmp = xAmp;

% time series of linear vegetation motion
t = (0:0.01:Tw)';
vege.t = t;
nt = length(t);
vege.x = zeros(nz,nt);
vege.z = zeros(nz,nt);
for i = 1:nt
    vege.x(:,i).....
        = Hv/2*GammaZ.*(gammaS.*sin(omegaW*t(i))+gammaC.*cos(omegaW*t(i)));
    vege.z(:,i) = s;
end

% total horizontal force on vegetation
theta=0:2*pi/100:2*pi;
nTheta = length(theta);
FI = zeros(nTheta,1);
FD = zeros(nTheta,1);
for iTheta = 1:nTheta
    FI(iTheta) = trapz(z, rho_w*Ac*Hv/2*omegaW^2*GammaZ....
        .*((1+Cm+Cm.*gammaS)*sin(theta(iTheta))....
        + (Cm.*gammaC)*cos(theta(iTheta))));
    FD(iTheta) = trapz(z, 1/2*Cd*rho_w*b.*(Hv/2*omegaW*GammaZ).^2.....
        .*abs((1+gammaS)*cos(theta(iTheta)) - gammaC*sin(theta(iTheta)))....
        .*((1+gammaS)*cos(theta(iTheta)) - gammaC*sin(theta(iTheta))));
end   %% be careful, Hv is used to consider incanopy velocity reduction
Fxrms = rms(FI+FD);
vege.Fxrms = Fxrms;

% wave decay
canopy.kD = kD;
if vege.tip
    xAmpTip2L = xAmp(1)/d2;
else
    xAmpTip2L = xAmp(end)/d2;
end
if canopy.alphaM4kD < 0
    if vege.tip
        if xAmpTip2L >= 0.015
            canopy.alphaM4kD = 0.2 * (xAmpTip2L)^(-0.5);
            % formula (32) and Figure 10(a) in Zhu et al.(2022)
            if xAmpTip2L >=0.155
                warning(['bottom xAmpTip2L = ',num2str(xAmpTip2L,'%5.3f'),....
                    ' formula may not be applicable, but the modification may be fine']);
            end
        else
            canopy.alphaM4kD = 1;
            warning(['bottom xAmpTip2L = ',num2str(xAmpTip2L,'%5.3f'),....
                ' formula for alphaM is not applicable so kD is not modified']);
        end
    else
        if xAmpTip2L >= 0.039
            canopy.alphaM4kD = 0.92 * (xAmpTip2L)^(-0.216);
            % formula (33) and Figure 10(b) in Zhu et al.(2022)
            if xAmpTip2L >= 0.449
                warning(['top xAmpTip2L = ',num2str(xAmpTip2L,'%5.3f'),....
                    ' formula may not be applicable, but the modification may be fine']);
            end
        else
            canopy.alphaM4kD = 1;
            warning(['top xAmpTip2L = ',num2str(xAmpTip2L,'%5.3f'),....
                ' formula for alphaM is not applicable so kD is not modified']);
        end
    end
end
canopy.kDModified = kD * canopy.alphaM4kD;
canopy.x = 0:Lv/20:Lv;
canopy.H = H0./(1 + canopy.kD*H0*canopy.x);
canopy.HModified = H0./(1 + canopy.kDModified*H0*canopy.x);

% bulk drag coefficient CD
canopy.CD = canopy.kD./canopy.alpha_eps./(b*N*k/(9*pi)....
    .* (9*sinh(k.*(d2+d3))-9*sinh(k.*d3) + sinh(3*k.*(d2+d3)) - sinh(3*k.*d3))....
    ./(sinh(k.*h).*(2*k.*h+sinh(2*k.*h)))); 
canopy.CDModified = canopy.kDModified./canopy.alpha_eps./(b*N*k/(9*pi)....
    .* (9*sinh(k.*(d2+d3))-9*sinh(k.*d3) + sinh(3*k.*(d2+d3)) - sinh(3*k.*d3))....
    ./(sinh(k.*h).*(2*k.*h+sinh(2*k.*h))));

%  effective vege length ratio Le2L = le/l
if vege.tip
    func = @(x) canopy.kD/canopy.alpha_eps - Cd(end).*b*N*k/(9*pi)....
        .* (9*sinh(k.*(h - d1))-9*sinh(k.*(h - d1 -x))...
        + sinh(3*k.*(h - d1)) - sinh(3*k.*(h - d1 -x)))....
        ./(sinh(k.*h).*(2*k.*h+sinh(2*k.*h)));
    canopy.Le2L = fzero(func, [0 d2])/d2;
    
    func = @(x) canopy.kDModified/canopy.alpha_eps - Cd(end).*b*N*k/(9*pi)....
        .* (9*sinh(k.*(h - d1))-9*sinh(k.*(h - d1 -x))...
        + sinh(3*k.*(h - d1)) - sinh(3*k.*(h - d1 -x)))....
        ./(sinh(k.*h).*(2*k.*h+sinh(2*k.*h)));
    canopy.Le2LModified = fzero(func, [0 d2])/d2;
else
    func = @(x) canopy.kD/canopy.alpha_eps - Cd(end).*b*N*k/(9*pi)....
        .* (9*sinh(k.*(x+d3))-9*sinh(k.*d3) + sinh(3*k.*(x+d3)) - sinh(3*k.*d3))....
        ./(sinh(k.*h).*(2*k.*h+sinh(2*k.*h)));
    canopy.Le2L = fzero(func, [0 d2])/d2;
    
    func = @(x) canopy.kDModified/canopy.alpha_eps - Cd(end).*b*N*k/(9*pi)....
        .* (9*sinh(k.*(x+d3))-9*sinh(k.*d3) + sinh(3*k.*(x+d3)) - sinh(3*k.*d3))....
        ./(sinh(k.*h).*(2*k.*h+sinh(2*k.*h)));
    canopy.Le2LModified = fzero(func, [0 d2])/d2;
end

%% plot results
% figure;clf;hold on;
% plot(vege.x/vege.d2, vege.z/vege.d2,'color',[0 0 0 0.25])
% axis equal    
% axis([-1 1 0 1])
% xlabel('Horizontal x/l')
% ylabel('Vertical z/l')
% title([{'Linear vegetation motion'},....
%     {'(see nonlinear vegetation motion in Zhu et al. 2020b and 2021)'}])
% 
% figure; clf; hold on;
% plot(canopy.x/canopy.Lv, canopy.H/H0,'b--')   
% plot(canopy.x/canopy.Lv, canopy.HModified/H0,'r-')  
% axis([0 1 0 1])
% xlabel('Distance x/L_v')
% ylabel('Wave height H/H_0')    
% title('Wave height decay over flexible vegetation canopy')
% legend('wave decay with linear vegetation motion',....
%     'modified wave decay by considering nonlinerity',....
%     'location','southwest')