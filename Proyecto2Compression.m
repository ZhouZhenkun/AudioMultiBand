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

% filtPlot = fvtool(genFilter1,1,genFilter2,1, genFilter3,1, genFilter4,1);
% legend(filtPlot,'Filter1','Filter2','Filter3','Filter4');
% filtPlot.Fs = Fs;

s1 = filter(genFilter1,1,inSignal);
s2 = filter(genFilter2,1,inSignal);
s3 = filter(genFilter3,1,inSignal);
s4 = filter(genFilter4,1,inSignal);
% 
% figure();
% t = (1:length(inSignal))/Fs;
% subplot(2,1,1)
% plot(t,inSignal)
% legend('Original Data')
% title('Separacion Bandas en Tiempo')
% subplot(2,1,2)
% plot(t,s1,t,s2,t,s3,t,s4)
% legend('Filter1','Filter2','Filter3','Filter4')
% 
% figure();
% xft = abs(fft(inSignal));
% 
% fx = (1:length(xft))*Fs;
% subplot(2,1,1)
% stem(fx,xft)
% legend('Original Data')
% title('Separacion Bandas en Frecuencia')
% subplot(2,1,2)
% sFilters = [abs(fft(s1)),abs(fft(s2)),abs(fft(s3)),abs(fft(s4))];
% stem(fx,sFilters)
% legend('Filter1','Filter2','Filter3','Filter4')


% Diezmado x4
sd1 = downsample(s1,numBands);
sd2 = downsample(s2,numBands);
sd3 = downsample(s3,numBands);
sd4 = downsample(s4,numBands);

% Apply reflection to s2 and s4 (odd) to invert their spectrum
% sd2

% 
% fd = (1:length(sd1)).*(Fs/numBands);
% fdn = (1:length(s1))/(Fs/numBands);
% 
% figure();
% subplot(2,1,1)
% stem(fdn,abs(fft(s1)))
% legend('Original Data')
% title('Separacion Bandas Downsample')
% subplot(2,1,2)
% stem(fd,abs(fft(sd1)))
% legend('Filter1')
% 
% figure();
% subplot(2,1,1)
% stem(fdn,abs(fft(s2)))
% legend('Original Data')
% title('Separacion Bandas Downsample')
% subplot(2,1,2)
% stem(fd,abs(fft(sd2)))
% legend('Filter2')
% 
% figure();
% subplot(2,1,1)
% stem(fdn,abs(fft(s3)))
% legend('Original Data')
% title('Separacion Bandas Downsample')
% subplot(2,1,2)
% stem(fd,abs(fft(sd3)))
% legend('Filter3')
% 
% figure();
% subplot(2,1,1)
% stem(fdn,abs(fft(s4)))
% legend('Original Data')
% title('Separacion Bandas Downsample')
% subplot(2,1,2)
% stem(fd,abs(fft(sd4)))
% legend('Filter4')

%% Compand

Mu = 255;
c1 = compand(sd1,Mu,max(sd1),'mu/compressor');
c2 = compand(sd2,Mu,max(sd2),'mu/compressor');
c3 = compand(sd3,Mu,max(sd3),'mu/compressor');
c4 = compand(sd4,Mu,max(sd4),'mu/compressor');

% Change this for scenarios
nCode1 =16;
nCode2 =1;
nCode3 =1;
nCode4 =1;

cCode1 = strcat('ubit',int2str(nCode1));
[partition,codebook] = lloyds(c1,2^nCode1);
[~,quants1] = quantiz(c1,partition,codebook);

cCode2 = strcat('ubit',int2str(nCode2));
[partition,codebook] = lloyds(c2,2^nCode2);
[~,quants2] = quantiz(c2,partition,codebook);

cCode3 = strcat('ubit',int2str(nCode3));
[partition,codebook] = lloyds(c3,2^nCode3);
[~,quants3] = quantiz(c3,partition,codebook);

cCode4 = strcat('ubit',int2str(nCode4));
[partition,codebook] = lloyds(c4,2^nCode4);
[~,quants4] = quantiz(c4,partition,codebook);


%%
nq1 = typecast(quants1,'uint64');
nq2 = typecast(quants2,'uint64');
nq3 = typecast(quants3,'uint64');
nq4 = typecast(quants4,'uint64');

shift1 = 64 - nCode1;
shift2 = 64 - nCode2;
shift3 = 64 - nCode3;
shift4 = 64 - nCode4;

nq1s = bitshift(nq1,-shift1);
nq2s = bitshift(nq2,-shift2);
nq3s = bitshift(nq3,-shift3);
nq4s = bitshift(nq4,-shift4);

%%

outBin = nq1s + nq2s + nq3s + nq4s;
fileID = fopen(outFileName,'w');
fwrite(fileID,nq1s, cCode1);
%fwrite(fileID,nq2s, cCode2);
%fwrite(fileID,nq3s, cCode3);
%fwrite(fileID,nq4s, cCode4);
fclose(fileID);

inFileStatus = dir(inputFileName);
outFileStatus = dir(outFileName);
fComp = inFileStatus.bytes/outFileStatus.bytes;
disp(inFileStatus.bytes)
disp(outFileStatus.bytes)
disp(fComp)


%% 
fileID = fopen(outFileName,'r');
dataIn1 = fread(fileID,length(nq1s),cCode1);
% dataIn2 = fread(fileID,length(nq2s),cCode2);
% dataIn3 = fread(fileID,length(nq3s),cCode3);
% dataIn4 = fread(fileID,length(nq4s),cCode4);
fclose(fileID);

dataInR1 = cast(dataIn1, 'uint64');
% dataInR2 = cast(dataIn2, 'uint64');
% dataInR3 = cast(dataIn3, 'uint64');
% dataInR4 = cast(dataIn4, 'uint64');
nquants1n = bitshift(dataInR1,shift1);
% nquants2n = bitshift(dataInR2,shift2);
% nquants3n = bitshift(dataInR3,shift3);
% nquants4n = bitshift(dataInR4,shift4);
quants1n = typecast(nquants1n, 'double');
% quants2n = typecast(nquants2n, 'double');
% quants3n = typecast(nquants3n, 'double');
% quants4n = typecast(nquants4n, 'double');


%%


newC1 = quants1n;
% newC2 = quants2n;
% newC3 = quants3n;
% newC4 = quants4n;

uc1 = compand(newC1,Mu,max(newC1),'mu/expander');
% uc2 = compand(newC2,Mu,max(newC2),'mu/expander');
% uc3 = compand(newC3,Mu,max(newC3),'mu/expander');
% uc4 = compand(newC4,Mu,max(newC4),'mu/expander');

newSd1 = upsample(uc1, numBands);
% newSd2 = upsample(uc2, numBands);
% newSd3 = upsample(uc3, numBands);
% newSd4 = upsample(uc4, numBands);

sn1 = filter(genFilter1,1,newSd1);
% sn2 = filter(genFilter2,1,newSd2);
% sn3 = filter(genFilter3,1,newSd3);
% sn4 = filter(genFilter4,1,newSd4);

%newRealSignal = (sn1 + sn2 + sn3 + sn4)*numBands;
newRealSignal = (sn1)*numBands;
audiowrite('outTest.wav',newRealSignal,Fs);

%% PESQ - Analysis

% name of executable file for PESQ calculation
binary = 'pesq2.exe';
% specify path to folder with reference and degraded audio files in it
pathaudio = '.';
% specify reference and degraded audio files
reference = inputFileName;
degraded = 'outTest.wav';

% compute NB-PESQ and WB-PESQ scores for wav-files
nb = pesq2_mtlb( reference, degraded, 8000, 'nb', binary, pathaudio );

% display results to screen
fprintf('====================\n'); 
disp('Compute NarrowBand-PESQ scores for wav-files:');
disp('(Top score = 4.5)');

MOS_percent = nb(1)*100/4.5;
MOS_LQO_percent = nb(2)*100/4.5;

% Mean Opinion Score
fprintf( 'NB PESQ MOS = %5.3f (%5.3f%%) \n', nb(1), MOS_percent);
% Mean Opinion Score - Liscening Quality Objective
fprintf( 'NB MOS LQO  = %5.3f (%5.3f%%) \n', nb(2), MOS_LQO_percent);
