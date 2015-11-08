obj.SampleRate = 44100; obj.PassbandRipple = 0.1; obj.StopbandAttenuation = 80;
obj.Bandwidth = 6000; obj.CenterFrequency = 4000;
Fs = obj.SampleRate;
Ap = obj.PassbandRipple;
Ast = obj.StopbandAttenuation;
Fp = obj.Bandwidth/2;
lpf = dsp.LowpassFilter('PassbandFrequency',Fp,...
'StopbandFrequency',Fp+1000,...
'PassbandRipple',Ap,...
'StopbandAttenuation',Ast,...
'SampleRate',Fs);
FIR = getFilter(lpf);
b = FIR.Numerator;
wc = obj.CenterFrequency*2/Fs;
L = size(b,2);
b = 2*cos(pi*wc*(-(L-1)/2:(L-1)/2)).*b;
FIR.Numerator = b;
obj.FilterObj = FIR;
fvtool(FIR,'Fs',Fs)



