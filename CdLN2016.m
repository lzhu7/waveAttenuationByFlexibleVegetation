function [Cd] = CdLN2016(Um, Tw, b)
% Drag coefficient of rectangular plates under waves calculated from the
% formulas in Luhar and Nepf (2016) and Luhar (2012,thesis) fitted from the
% experimental data from Keulegan and Carpenter (1958) and Sarpkaya and
% O’Keefe (1996).

%% input variables
% Um: [m/s] amplitude of relative velocity;
% Tw: [s] wave period;
% b: [m] blade width projected to the direction of wave propagation

%% output variables
% Cd: [-] added mass coefficient

%%
KC = Um*Tw/b;
Cd = max(10*KC.^(-1/3),1.95);