function [Cm] = CmLN2016(Um, Tw, b)
% Added mass coefficient of rectangular plates under waves calculated from
% the formulas in Luhar and Nepf (2016) and Luhar (2012,thesis) fitted from
% the experimental data from Keulegan and Carpenter (1958) and Sarpkaya and
% O’Keefe (1996).

%% input variables
% Um: [m/s] amplitude of relative velocity;
% Tw: [s] wave period;
% b: [m] blade width projected to the direction of wave propagation

%% output variables
% Cm: [-] added mass coefficient

%%
KC = Um*Tw/b;
nz = length(Um);
Cm = zeros(nz,1);
for i=1:nz
    if KC(i)<20
        CM1 = 1+0.35*KC(i)^(2/3);
    else
        CM1 = 1+0.15*KC(i)^(2/3);
    end
    CM2 = 1+(KC(i)-18)^2/49;
    Cm(i) = min(CM1,CM2);
end
