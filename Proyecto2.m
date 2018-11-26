close all;
clear all;
clc;
%% Initial Parameters
inputFileName = 'test2.wav';
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
fdn = (1:length(s1))/(Fs/numBands);

figure();
subplot(2,1,1)
stem(fdn,abs(fft(s1)))
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd1)))
legend('Filter1')

figure();
subplot(2,1,1)
stem(fdn,abs(fft(s2)))
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd2)))
legend('Filter2')

figure();
subplot(2,1,1)
stem(fdn,abs(fft(s3)))
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd3)))
legend('Filter3')

figure();
subplot(2,1,1)
stem(fdn,abs(fft(s4)))
legend('Original Data')
title('Separacion Bandas Downsample')
subplot(2,1,2)
stem(fd,abs(fft(sd4)))
legend('Filter4')

%% Compand

Mu = 255;
c1 = compand(sd1,Mu,max(sd1),'mu/compressor');
c2 = compand(sd2,Mu,max(sd2),'mu/compressor');
c3 = compand(sd3,Mu,max(sd3),'mu/compressor');
c4 = compand(sd4,Mu,max(sd4),'mu/compressor');




binComp1 = typecast(c1, 'uint16');
binComp2 = typecast(c2, 'uint16');
binComp3 = typecast(c3, 'uint16');
binComp4 = typecast(c4, 'uint16');

%[~,~,distor] = quantiz(sd1,partition,codebook);
nCode =15;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants1] = quantiz(binComp1,partition,codebook);
binData1= dec2bin(quants1,nCode);

nCode =10;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants2] = quantiz(binComp2,partition,codebook);
binData2= dec2bin(quants2,nCode);

nCode =8;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants3] = quantiz(binComp3,partition,codebook);
binData3= dec2bin(quants3,nCode);

nCode =3;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants4] = quantiz(binComp4,partition,codebook);
binData4= dec2bin(quants4,nCode);

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
outBin = [binData1, binData2, binData3, binData4];
fileID = fopen(outFileName,'w');
fwrite(fileID,outBin, 'ubit40');
fclose(fileID);

inFileStatus = dir(inputFileName);
outFileStatus = dir(outFileName);
fComp = inFileStatus.bytes/outFileStatus.bytes;
disp(inFileStatus.bytes)
disp(outFileStatus.bytes)
disp(fComp)
