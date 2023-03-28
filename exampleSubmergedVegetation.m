clear; clc; close all;
% This is an example for submerged, i.e., case H1 in Luhar et al. (2017)
% and Figure 5(d) in Zhu et al. (2022)

%% seagrass blade (flexible)
% input variables and units are defined in the function
% regularWaveDecayFlexibleVegetation (type 'help
% regularWaveDecayFlexibleVegetation' to find detail information)
%
% wave parameters
    wave.h = 0.16;      % [m] water depth
    wave.rho = 1000;   % [k/m^3] water mass density
    wave.Tw = 1.4286;       % [s] wave period
    wave.H0 = 0.028;     % [m] wave height
% blade parameters
    blade.rho = 920;    % [kg/m^3] blade mass density
    blade.EI = 7.5000e-08;  % [Nm^2] blade flexural rigidity                       
    blade.b = 0.0030;     % [m] blade width
    blade.Ac = 3.0000e-07;    % [m^2] area of blade cross section
    blade.d2 = 0.1300;    % [m] blade length
    blade.d3 = 0.0100;       % [m] gap below the canopy
    blade.tip = 0;      % 0 --> fixed at bottom; 1 --> fixed at top
% blade canopy parameters
    bladeCanopy.Lv = 5;       % [m] canopy length in wave direction
    bladeCanopy.N = 7200;      % [baldes/m^2] canopy/plant density
    bladeCanopy.alpha_eps = 1;  % sheltering factor for blades
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
                                % 
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

%% sheath (rigid)
% input variables and units are defined in the function
% regularWaveDecayRigidVegetation (type 'help
% regularWaveDecayRigidVegetation' to find detail information)
%
% wave parameters (use the same values from flexbile blade)
% sheath parameters
    sheath.b = 0.0078;     % [m] stem width
    sheath.Ac = 4.7784e-05;    % [m^2] area of stem cross section
    sheath.d2 = 0.0100;    % [m] stem length
    sheath.d3 = 0;       % [m] gap below the stem canopy
    sheath.tip = 0;      % 0 --> fixed at bottom; 1 --> fixed at top
% sheath canopy parameters
    sheathCanopy.Lv = 5;       % [m] canopy length in wave direction
    sheathCanopy.N = 1200;      % [baldes/m^2] canopy/plant density
    sheathCanopy.alpha_w = -1; % default value -1
                      
% sheath hydrodynamic coefficients
    sheathHydroCoeff.Cd = -1; % default value -1
    sheathHydroCoeff.Cm = -1; % default value -1
    
% calculate wave attenuation by rigid part
[sheath, sheathCanopy, sheathHydroCoeff]....
= regularWaveDecayRigidVegetation(wave, sheath, sheathCanopy, sheathHydroCoeff);
    
%% wave attenuation by the whole kelp canopy
seagrassMeadow.kD = bladeCanopy.kD + sheathCanopy.kD;
seagrassMeadow.kDModified = bladeCanopy.kDModified + sheathCanopy.kD;
seagrassMeadow.x = bladeCanopy.x;
seagrassMeadow.Lv = bladeCanopy.Lv;
seagrassMeadow.H = wave.H0./(1 + seagrassMeadow.kD*wave.H0*seagrassMeadow.x);
seagrassMeadow.HModified = wave.H0./(1 + seagrassMeadow.kDModified*wave.H0*seagrassMeadow.x);

%% plot results
% blade motion
figure;clf;hold on;
plot(blade.x/blade.d2, blade.z/blade.d2,'color',[0 0 0 0.25])
axis equal    
axis([-1 1 0 1])
xlabel('Horizontal x/l')
ylabel('Vertical z/l')
title([{'Linear motion of a submerged seagrass blade'},....
    {'(see nonlinear motion in Zhu et al. 2020, JGR Oceans)'}])

% wave decay
% figureElsevier11
figure;clf; hold on
data = load('seagrassH1.mat');
plot(data.x/seagrassMeadow.Lv, data.H2H0,'k.') 
plot(seagrassMeadow.x/seagrassMeadow.Lv, seagrassMeadow.H/wave.H0,'b--')   
plot(seagrassMeadow.x/seagrassMeadow.Lv, seagrassMeadow.HModified/wave.H0,'r')   
axis([0 1 0 1.1])
xlabel('Distance x/L_v')
ylabel('Wave height H/H_0')    
title('Wave height decay over submerged seagrass meadow')
legend('measured data for case H1 in Luhar et al. (2017)','wave decay with linear blade motion',....
    'modified wave decay by considering nonlinerity',....
    'location','southwest')
    
