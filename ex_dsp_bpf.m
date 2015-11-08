sw = dsp.SineWave('SampleRate',44.1e3,'Frequency',[100,10000],...
                'SamplesPerFrame',1024);
bpf = dsp.BandpassFilter;
ap = dsp.ArrayPlot;
tic,
while toc < 10
    in = sum(step(sw),2);
    out = step(bpf,in);
    step(ap,[in out]);
end