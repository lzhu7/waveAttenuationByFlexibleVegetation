%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2022%%%%%%%%%%%%%%%%%%%%%%
%---------------------------version 2022.04--------------------------------
%
% ./
% 
% This is a package of MATLAB transcripts for the analytical models for
% regular/monochromatic wave attenuation by submerged and suspended
% flexible and rigid vegetation and kelp (Zhu and Zou, 2017; Zhu et al.,
% 2022), including two major functions,
% regularWaveDecayFlexibleVegetation.m and
% regularWaveDecayRigidVegetation.m, respectively, for wave attenuation by
% flexible and rigid vegetation/kelp. Type 'help' and the name of the
% function to see details.
% 
% Two EXAMPLES, exampleSubmergedVegetation.m and exampleSuspendedKelp.m,
% are given for submerged vegetation and suspended kelp, respectively.
%
% Recommended Citation:
% [1] Zhu, L., & Zou, Q. (2017). Three-layer analytical solution for wave
% attenuation by suspended and nonsuspended vegetation canopy. Coast. Eng.
% Proc. 1 (35), 27. https://doi.org/ 10.9753/icce.v35.waves.27.
%
% [2] Zhu, L., Huguenard, K., Fredriksson, D. W., & Lei, J. (2022). Wave
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
% [1] Zhu, L., & Zou, Q. (2017). Three-layer analytical solution for wave
% attenuation by suspended and nonsuspended vegetation canopy. Coast. Eng.
% Proc. 1 (35), 27. https://doi.org/ 10.9753/icce.v35.waves.27.
%
% [2] Zhu, L. (2020). Wave Attenuation Capacity of Suspended Aquaculture
% Structures with Sugar Kelp and Mussels. PhD Thesis. The University of
% Maine. https://digitalcommons.library.umaine.edu/etd/3222
%
% [3] Zhu, L., Huguenard, K., Zou, Q., Fredriksson, D.W., & Xie, D.
% (2020a). Aquaculture farms as nature-based coastal protection: random
% wave attenuation by suspended and submerged canopies. Coast. Eng. 160,
% 103737 https://doi.org/10.1016/j.coastaleng.2020.103737
%
% [4] Zhu, L., Zou, Q., Huguenard, K., & Fredriksson, D. W. (2020b).
% Mechanisms for the Asymmetric Motion of Submerged Aquatic Vegetation in
% Waves: A Consistent-Mass Cable Model. JGR: Oceans, 125(2).
% https://doi.org/10.1029/2019JC015517
%
% [5] Zhu, L., Lei, J., Huguenard, K., & Fredriksson, D.W. (2021). Wave
% attenuation by suspended canopies with cultivated kelp (Saccharina
% latissima). Coast. Eng., 103947
% https://doi.org/10.1016/j.coastaleng.2021.103947
%
% [6] Zhu, L., Huguenard, K., Fredriksson, D. W., & Lei, J. (2022). Wave
% attenuation by flexible vegetation (and suspended kelp) with blade
% motion: Analytical solutions. Advances in Water Resources, 162, 104148.
% https://doi.org/10.1016/j.advwatres.2022.104148
%
% Contact longhuan.zhu@maine.edu in case of any questions.
%
%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2022%%%%%%%%%%%%%%%%%%%%%%
%
