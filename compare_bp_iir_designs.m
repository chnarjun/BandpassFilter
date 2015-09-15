
%% Default constructor
bp = fdesign.bandpass;
hbut = design(bp,'butter');
hch1 = design(bp,'cheby1');
hch2 = design(bp,'cheby2');
hell = design(bp,'ellip');
fvtool(hbut,hch1,hch2,hell);
legend('Butter','cheby1','cheby2','Ellip');

%% Wide transition band
bp = fdesign.bandpass(0.01,0.45,0.55,0.99);
hbut = design(bp,'butter'); % Far away from specs
hch1 = design(bp,'cheby1');
hch2 = design(bp,'cheby2'); % Far away from specs
hell = design(bp,'ellip');
fvtool(hbut,hch1,hch2,hell);
legend('Butter','cheby1','cheby2','Ellip');

%% Narrow transition band
bp = fdesign.bandpass(0.01,0.49,0.51,0.99);
hbut = design(bp,'butter');  % Far away from specs
hch1 = design(bp,'cheby1');  
hch2 = design(bp,'cheby2'); % Far away from specs
hell = design(bp,'ellip');
fvtool(hbut,hch1,hch2,hell);
legend('Butter','cheby1','cheby2','Ellip');

%% Narrow transition band, very wide passband 
% Pretty much all designs met the specs
bp = fdesign.bandpass(0.01,0.02,0.98,0.99);
hbut = design(bp,'butter');  
hch1 = design(bp,'cheby1');  
hch2 = design(bp,'cheby2'); 
hell = design(bp,'ellip');
fvtool(hbut,hch1,hch2,hell);
legend('Butter','cheby1','cheby2','Ellip');


%% NOTES

% For IIR Bandpass design go with Cheby1 design as it almost matches
% specifications and cost of filter design is less for all the cases. Next
% best in terms of cost is elliptical. 


% Minimum order design that takes Fst1, Fp1, Fp2, Fst2,Ast1, Ap, Ast2