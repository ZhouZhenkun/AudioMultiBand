close all
%% Initial Parameters
[inSignal, Fs] = audioread('testTone.wav');
%[inSignal, Fs] = audioread('AudioDePrueba.wav');

numBands = 4;

filters = zeros(numBands);
%% Symmetric Sub-Band Decimation
rp = 1;
rs = 60;
devPass = (10^(rp/20)-1)/(10^(rp/20)+1);
devStop = 10^(-rs/20);

fc = (Fs/2)/numBands;
f = [fc fc*1.1];
a= [1, 0];
dev = [devPass, devStop ]; 
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter1 = firpm(n,fo,ao,w);

f = [fc*0.9, fc, 2*fc, 2*fc*1.1];
a= [0, 1, 0];
dev = [devStop, devPass, devStop];  
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter2 = firpm(n,fo,ao,w);

f = [2*fc*0.9, 2*fc, 2*fc + fc, (2*fc + fc)*1.1];
a= [0, 1, 0];
dev = [devStop, devPass, devStop];  
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter3 = firpm(n,fo,ao,w);


fc = 3*fc;
f = [0.9*fc, fc];
a= [0, 1];
dev = [devStop, devPass]; 
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter4 = firpm(n,fo,ao,w);

filtPlot = fvtool(genFilter1,1,genFilter2,1, genFilter3,1, genFilter4,1);
legend(filtPlot,'Filter1','Filter4','Filter3','Filter4');
filtPlot.Fs = Fs;

s1 = filter(genFilter1,1,inSignal);
s2 = filter(genFilter2,1,inSignal);
s3 = filter(genFilter3,1,inSignal);
s4 = filter(genFilter4,1,inSignal);

figure();
t = (1:length(inSignal))/Fs;
subplot(2,1,1)
plot(t,inSignal)
legend('Original Data')
title('Separacion Bandas')
subplot(2,1,2)
plot(t,s1,t,s2,t,s3,t,s4)
legend('Filter1','Filter2','Filter3','Filter4')

figure();
xft = abs(fft(inSignal));

fx = (1:length(xft))*Fs;
subplot(2,1,1)
stem(fx,xft)
legend('Original Data')
title('Separacion Bandas')
subplot(2,1,2)
sFilters = [abs(fft(s1)),abs(fft(s2)),abs(fft(s3)),abs(fft(s4))];
stem(fx,sFilters)
legend('Filter1','Filter2','Filter3','Filter4')
%% Compand

Mu = 255;