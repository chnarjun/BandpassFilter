
%% Default constructor
bp = fdesign.bandpass;
heqrip = design(bp,'equiripple');
hkaiser = design(bp,'kaiserwin');
fvtool(heqrip,hkaiser);
legend('Equiripple','Kaiser');

%% Wide transition band
bp = fdesign.bandpass(0.01,0.45,0.55,0.99);
heqrip = design(bp,'equiripple'); % Doesn't converge
hkaiser = design(bp,'kaiserwin');
fvtool(heqrip,hkaiser);
legend('Equiripple','Kaiser');

%% Narrow transition band
bp = fdesign.bandpass(0.01,0.49,0.51,0.99);
heqrip = design(bp,'equiripple');
hkaiser = design(bp,'kaiserwin');
fvtool(heqrip,hkaiser);
legend('Equiripple','Kaiser');

%% Narrow transition band, very wide passband 
bp = fdesign.bandpass(0.01,0.02,0.98,0.99);
heqrip = design(bp,'equiripple');
hkaiser = design(bp,'kaiserwin');
fvtool(heqrip,hkaiser);
legend('Equiripple','Kaiser');


%%  NOTES

% For FIR design of Bandpass filter, it is better to go with kaiser window
% as equiripple design doesn't converge for large transition widths. See
% case 2 above
