close all;
clear all;
clc;
%% Initial Parameters
inputFileName = 'test.wav';
outFileName = 'audio.bin';
[inSignal, Fs] = audioread(inputFileName, 'native');
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
sd2Mirror = downsample(s2,numBands);
sd3 = downsample(s3,numBands);
sd4Mirror = downsample(s4,numBands);

% Apply reflection to s2 and s4 (odd) to invert their spectrum
% sd2
for (i = 1:length(sd2Mirror))
  if (mod(i,2) == 0)
    sd2(i) = sd2Mirror(i);
  else 
    sd2(i) = -sd2Mirror(i);
  end
end
sd2 = transpose(sd2);

% sd4
for (i = 1:length(sd4Mirror))
  if (mod(i,2) == 0)
    sd4(i) = sd4Mirror(i);
  else 
    sd4(i) = -sd4Mirror(i);
  end
end
sd4 = transpose(sd4);

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

binComp1 = typecast(round(c1), 'uint16');
binComp2 = typecast(round(c2), 'uint16');
binComp3 = typecast(round(c3), 'uint16');
binComp4 = typecast(round(c4), 'uint16');

% binComp1 = round(c1);
% binComp2 = round(c2);
% binComp3 = round(c3);
% binComp4 = round(c4);

%[~,~,distor] = quantiz(sd1,partition,codebook);
nCode =15;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants1] = quantiz(binComp1,partition,codebook);
%binData1= de2bi(quants1,nCode+1);

nCode =9;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants2] = quantiz(binComp2,partition,codebook);
%binData2= de2bi(quants2,nCode+1);

nCode =7;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants3] = quantiz(binComp3,partition,codebook);
%binData3= de2bi(quants3,nCode+1);

nCode =5;
partition = 0:2^nCode-1;
codebook =  0:2^nCode;
[~,quants4] = quantiz(binComp4,partition,codebook);
%binData4= de2bi(quants4,nCode+1);

%%
%outBin = [binData1, binData2, binData3, binData4];
nq1 = cast(quants1, 'uint16');
nq2 = cast(quants2, 'uint16');
nq3 = cast(quants3, 'uint16');
nq4 = cast(quants4, 'uint16');

outBin = [nq1, nq2, nq3, nq4];
fileID = fopen(outFileName,'w');
fwrite(fileID,outBin);
% fwrite(fileID,binData2, 'ubit14');
% fwrite(fileID,binData3, 'ubit12');
% fwrite(fileID,binData4, 'ubit8');
fclose(fileID);

inFileStatus = dir(inputFileName);
outFileStatus = dir(outFileName);
fComp = inFileStatus.bytes/outFileStatus.bytes;
disp(inFileStatus.bytes)
disp(outFileStatus.bytes)
disp(fComp)


%% 
fileID = fopen(outFileName,'r');
dataIn = fread(fileID,outFileStatus.bytes,'uint16');
fclose(fileID);

%%

newData1 = nq1;
newC1C = typecast(newData1, 'double');
newC1 = round(newC1C);

newData2 = nq2;
newC2C = cast(newData2, 'uint16');
newC2 = typecast(newData2, 'double');

newData3 = nq3;
newC3C = cast(newData3, 'uint16');
newC3 = typecast(newData3, 'double');

newData4 = nq4;
newC4C = cast(newData4, 'uint16');
newC4 = typecast(newData4, 'double');

uc1 = compand(newC1,Mu,max(newC1),'mu/expander');
uc2 = compand(newC2,Mu,max(newC2),'mu/expander');
uc3 = compand(newC3,Mu,max(newC3),'mu/expander');
uc4 = compand(newC4,Mu,max(newC4),'mu/expander');

newSd1 = upsample(uc1, numBands);
newSd2 = upsample(uc2, numBands);
newSd3 = upsample(uc3, numBands);
newSd4 = upsample(uc4, numBands);

sn1 = filter(genFilter1,1,newSd1);
sn2 = filter(genFilter2,1,newSd2);
sn3 = filter(genFilter3,1,newSd3);
sn4 = filter(genFilter4,1,newSd4);

newRealSignal = (sn1 + sn2 + sn3 + sn4)*numBands;

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







