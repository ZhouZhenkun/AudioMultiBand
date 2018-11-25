close all;
clear all;
clc;
%% Initial Parameters
inputFileName = 'test.wav';
outFileName = 'audio.bin';
[inSignal, Fs] = audioread(inputFileName);
%[inSignal, Fs] = audioread('AudioDePrueba.wav');

numBands = 4;
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
legend(filtPlot,'Filter1','Filter2','Filter3','Filter4');
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
title('Separacion Bandas en Tiempo')
subplot(2,1,2)
plot(t,s1,t,s2,t,s3,t,s4)
legend('Filter1','Filter2','Filter3','Filter4')

figure();
xft = abs(fft(inSignal));

fx = (1:length(xft))*Fs;
subplot(2,1,1)
stem(fx,xft)
legend('Original Data')
title('Separacion Bandas en Frecuencia')
subplot(2,1,2)
sFilters = [abs(fft(s1)),abs(fft(s2)),abs(fft(s3)),abs(fft(s4))];
stem(fx,sFilters)
legend('Filter1','Filter2','Filter3','Filter4')


sd1 = downsample(s1,numBands);
sd2 = downsample(s2,numBands);
sd3 = downsample(s3,numBands);
sd4 = downsample(s4,numBands);


fd = (1:length(sd1)).*(Fs/numBands);
td = (1:length(sd1))/(Fs/numBands);

figure();
subplot(2,1,1)
stem(td,sd1)
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd1)))
legend('Filter1')

figure();
subplot(2,1,1)
stem(td,sd2)
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd2)))
legend('Filter2')

figure();
subplot(2,1,1)
stem(td,sd3)
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd3)))
legend('Filter3')

figure();
subplot(2,1,1)
stem(td,sd4)
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd4)))
legend('Filter4')

%% Compand

Mu = 255;
c1 = compand(sd1,Mu,max(sd1),'mu/compressor');
c2 = compand(sd1,Mu,max(sd2),'mu/compressor');
c3 = compand(sd1,Mu,max(sd3),'mu/compressor');
c4 = compand(sd1,Mu,max(sd4),'mu/compressor');


partition = 0:2^6-1;
codebook = 0:2^6;
%[~,~,distor] = quantiz(sd1,partition,codebook);
[~,quants1] = quantiz(c1,partition,codebook);
[~,quants2] = quantiz(c2,partition,codebook);
[~,quants3] = quantiz(c3,partition,codebook);
[~,quants4] = quantiz(c4,partition,codebook);


% newsig = compand(quants1,Mu,max(quants1),'mu/expander');
% 
% newsig2 = reshape(newsig, [length(newsig),1]);
% distor2 = sum((newsig2-sd1).^2)/length(sd1);
% 
% figure();
% subplot(2,1,1)
% plot(sd1)
% legend('Original Data')
% title('Comparacion')
% subplot(2,1,2)
% plot(c1)
% legend('Companded')

%%
outBin = [quants1, quants2, quants3, quants4];
fileID = fopen(outFileName,'w');
fwrite(fileID,outBin);
fclose(fileID);

inFileStatus = dir(inputFileName);
outFileStatus = dir(outFileName);
fComp = inFileStatus.bytes./outFileStatus.bytes;