bp = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',1000);
hfirls = design(bp,'firls');
heqrip = design(bp,'equiripple');
fvtool(hfirls,heqrip);

% FIRCLS
bp = fdesign.bandpass('N,Fc1,Fc2,Ast1,Ap,Ast2',1000);
hfircls = design(bp,'fircls');
fvtool(hfircls);

% ELLIP 
bp = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2' ,126);
hellip = design(bp,'ellip');
fvtool(hellip);