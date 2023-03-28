clear; clc; close all;
% This is an example for suspended kelp, i.e., case 5 in Zhu et al. (2021)
% and Figure 5(a) in Zhu et al. (2022)

%% kelp blade (flexible)
% input variables and units are defined in the function
% regularWaveDecayFlexibleVegetation (type 'help
% regularWaveDecayFlexibleVegetation' to find detail information)
%
% wave parameters
    wave.h = 0.4;      % [m] water depth
    wave.rho = 1000;   % [k/m^3] water mass density
    wave.Tw = 1;       % [s] wave period
    wave.H0 = 0.039;     % [m] wave height
% blade parameters
    blade.rho = 1200;    % [kg/m^3] blade mass density
    blade.EI = 1.6176e-9;  % [Nm^2] blade flexural rigidity                       
    blade.b = 0.0095;     % [m] blade width
    blade.Ac = 9.5250e-07;    % [m^2] area of blade cross section
    blade.d2 = 0.0966;    % [m] blade length
    blade.d3 = 0.1884;       % [m] gap below the canopy
    blade.tip = 1;      % 0 --> fixed at bottom; 1 --> fixed at top
% blade canopy parameters
    bladeCanopy.Lv = 3.8;       % [m] canopy length in wave direction
    bladeCanopy.N = 5.2632e+03;      % [baldes/m^2] canopy/plant density
    bladeCanopy.alpha_eps = 0.6301;  % sheltering factor for blades
    bladeCanopy.alpha_w = -1; % default value -1: calculate alpha_w the
%                            % formula in Lowe et al. (2005)
%                            % set a positive value if it is known
%                            % e.g., 1 --> no velocity reduction in canopy
    bladeCanopy.alphaM4kD = -1; % default value -1: calculate alpha_M using
%                            % (32) and (33) in Zhu et al. (2022) for
%                            % submerged vegetation and suspended
%                            % vegetation, respectively.
%                            % set a positive value if it is known
%                            % e.g., 1 --> no modification                       
% blade hydrodynamic coefficients
    bladeHydroCoeff.Cd = -1; % default value -1: calcualte Cd using the
%                            % formula in Luhar and Nepf (2016).
%                            % set a positive value if it is known
%                            % e.g., 1.95 --> drag coefficient is 1.95
    bladeHydroCoeff.Cm = -1; % default value -1: calcualte Cm using the
%                            % formula in Luhar and Nepf (2016).
%                            % set a positive value if it is known
%                            % e.g., 1 --> added mass coefficient is 1
    
% calculate blade motion and wave attenuation by flexible blades
[blade, bladeCanopy, bladeHydroCoeff] ...
    = regularWaveDecayFlexibleVegetation(wave,blade,bladeCanopy,bladeHydroCoeff);

%% kelp stipe (rigid)
% input variables and units are defined in the function
% regularWaveDecayRigidVegetation (type 'help
% regularWaveDecayRigidVegetation' to find detail information)
%
% wave parameters (input already)
% stipe parameters
    stipe.b = 0.0095;     % [m] stem width
    stipe.Ac = 9.5250e-06;    % [m^2] area of stem cross section
    stipe.d2 = 0.0050;    % [m] stem length
    stipe.d3 = 0.2850;       % [m] gap below the stem canopy
    stipe.tip = 1;      % 0 --> fixed at bottom; 1 --> fixed at top
% stipe canopy parameters
    stipeCanopy.Lv = 3.8;       % [m] canopy length in wave direction
    stipeCanopy.N = 526.3158;      % [baldes/m^2] canopy/plant density
    stipeCanopy.alpha_w = -1; % default value -1
                      
% stipe hydrodynamic coefficients
    stipeHydroCoeff.Cd = -1; % default value -1
    stipeHydroCoeff.Cm = -1; % default value -1
    
% calculate wave attenuation by rigid part
[stipe, stipeCanopy, stipeHydroCoeff]....
= regularWaveDecayRigidVegetation(wave, stipe, stipeCanopy, stipeHydroCoeff);
    
%% wave attenuation by the whole kelp canopy
kelpCanopy.kD = bladeCanopy.kD + stipeCanopy.kD;
kelpCanopy.kDModified = bladeCanopy.kDModified + stipeCanopy.kD;
kelpCanopy.x = bladeCanopy.x;
kelpCanopy.Lv = bladeCanopy.Lv;
kelpCanopy.H = wave.H0./(1 + kelpCanopy.kD*wave.H0*kelpCanopy.x);
kelpCanopy.HModified = wave.H0./(1 + kelpCanopy.kDModified*wave.H0*kelpCanopy.x);

%% plot results
% blade motion
figure;clf;hold on;
plot(blade.x/blade.d2, blade.z/blade.d2,'color',[0 0 0 0.25])
axis equal    
axis([-1 1 0 1])
xlabel('Horizontal x/l')
ylabel('Vertical z/l')
title([{'Linear blade motion'},....
    {'(see nonlinear motion in Zhu et al. 2021, Coastal Engineering)'}])

% wave decay
% figureElsevier11
figure;clf; hold on
data = load('kelpCase5.mat');
plot(data.x/kelpCanopy.Lv, data.H2H0,'k.')  
plot(kelpCanopy.x/kelpCanopy.Lv, kelpCanopy.H/wave.H0,'b--')   
plot(kelpCanopy.x/kelpCanopy.Lv, kelpCanopy.HModified/wave.H0,'r')   
axis([0 1 0 1.1])
xlabel('Distance x/L_v')
ylabel('Wave height H/H_0')    
title('Wave height decay over suspended kelp canopy')
legend('measured data for Case 5 in Zhu et al. (2021)','wave decay with linear blade motion',....
    'modified wave decay by considering nonlinerity',....
    'location','southwest')