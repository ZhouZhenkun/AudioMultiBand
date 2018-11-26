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

% 1dB = 0.5dB during decimation + 0.5dB during interpolation
rp = 1;                                                         % 1dB? no deberia de ser 0.5dB?
% 60dB attenuation
rs = 60;
devPass = (10^(rp/20)-1)/(10^(rp/20)+1);
devStop = 10^(-rs/20);

% Filter1 - LowPass fc = Fs/8
fc = (Fs/2)/numBands;
f = [fc fc*1.1];
a= [1, 0];
dev = [devPass, devStop ]; 
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter1 = firpm(n,fo,ao,w);

% Filter2 - PassBand fc1 = Fs/8, fc2 = Fs/4
f = [fc*0.9, fc, 2*fc, 2*fc*1.1];
a= [0, 1, 0];
dev = [devStop, devPass, devStop];  
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter2 = firpm(n,fo,ao,w);

% Filter3 - PassBand fc1 = Fs/4, fc2 = 3*Fs/8
f = [2*fc*0.9, 2*fc, 2*fc + fc, (2*fc + fc)*1.1];
a= [0, 1, 0];
dev = [devStop, devPass, devStop];  
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter3 = firpm(n,fo,ao,w);

% Filter4 - HighPass fc = 3*Fs/8
fc = 3*fc;
f = [0.9*fc, fc];
a= [0, 1];
dev = [devStop, devPass]; 
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
genFilter4 = firpm(n,fo,ao,w);

% Display all 4 filters in the Filter Visualization Tool
filtPlot = fvtool(genFilter1,1,genFilter2,1, genFilter3,1, genFilter4,1);
legend(filtPlot,'Filter1 (LowPass)','Filter2 (PassBand)','Filter3 (PassBand)','Filter4 (HighPass)');
filtPlot.Fs = Fs;

% Apply Filtering to input signal
s1 = filter(genFilter1,1,inSignal);
s2 = filter(genFilter2,1,inSignal);
s3 = filter(genFilter3,1,inSignal);
s4 = filter(genFilter4,1,inSignal);

% Plot results in time
figure();
t = (1:length(inSignal))/Fs;
subplot(2,1,1)
plot(t,inSignal)
ylim([-1 1])
legend('Original Data')
title('Separacion Bandas en Tiempo')
subplot(2,1,2)
plot(t,s1,t,s2,t,s3,t,s4)
ylim([-1 1])
legend(filtPlot,'Filter1 (LowPass)','Filter2 (PassBand1)','Filter3 (PassBand2)','Filter4 (HighPass)');

% Plot results in frequency
% figure();
% xft = abs(fft(inSignal));
% fx = (1:length(xft))*Fs;
% subplot(2,1,1)
% stem(fx,xft)
% legend('Original Data')
% title('Separacion Bandas en Frecuencia')
% subplot(2,1,2)
% sFilters = [abs(fft(s1)),abs(fft(s2)),abs(fft(s3)),abs(fft(s4))];
% stem(fx,sFilters)
% legend(filtPlot,'Filter1 (LowPass)','Filter2 (PassBand1)','Filter3 (PassBand2)','Filter4 (HighPass)');

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

% Filter undesired alias images using LowPass filter with fc = Fs/8
% (reusing Lowpass filter previously designed)
sdf1 = filter(genFilter1,1,sd1);
sdf2 = filter(genFilter1,1,sd2);
sdf3 = filter(genFilter1,1,sd3);                % Porque todavia veo info alrededor de 20kHz?
sdf4 = filter(genFilter1,1,sd4);

% Plot the downsampled results for all 4 segments
fd = (1:length(sd1)).*(Fs/numBands);
td = (1:length(sd1))/(Fs/numBands);

% First segment
figure();
subplot(2,2,1)
stem(td,sd1)
legend('Original Data')
title('Separacion Bandas Downsample - sd1')
subplot(2,2,2)
stem(fd,abs(fft(sd1)))
legend('Filter1 (LowPass)')

subplot(2,2,3)
stem(td,sdf1)
legend('Original Data')
title('Separacion Bandas Downsample - sd1 filtrada')
subplot(2,2,4)
stem(fd,abs(fft(sdf1)))
legend('Filter1 (LowPass)')

% 2nd segment
% figure();
% subplot(2,2,1)
% stem(td,sd2)
% legend('Original Data')
% title('Separacion Bandas Downsample - sd2')
% subplot(2,2,2)
% stem(fd,abs(fft(sd2)))
% legend('Filter2 (PassBand1)')
% 
% subplot(2,2,3)
% stem(td,sdf2)
% legend('Original Data')
% title('Separacion Bandas Downsample - sd2 filtrada')
% subplot(2,2,4)
% stem(fd,abs(fft(sdf2)))
% legend('Filter2 (PassBand1)')

% 3rd segment
% figure();
% subplot(2,2,1)
% stem(td,sd3)
% legend('Original Data')
% title('Separacion Bandas Downsample - sd3')
% subplot(2,2,2)
% stem(fd,abs(fft(sd3)))
% legend('Filter3 (PassBand2)')
% 
% subplot(2,2,3)
% stem(td,sdf3)
% legend('Original Data')
% title('Separacion Bandas Downsample - sd3 filtrada')
% subplot(2,2,4)
% stem(fd,abs(fft(sdf3)))
% legend('Filter3 (PassBand2)')

% 4rd segment
% figure();
% subplot(2,2,1)
% stem(td,sd4)
% legend('Original Data')
% title('Separacion Bandas Downsample - sd4')
% subplot(2,2,2)
% stem(fd,abs(fft(sd4)))
% legend('Filter4 (HighPass)')
% 
% subplot(2,2,3)
% stem(td,sdf4)
% legend('Original Data')
% title('Separacion Bandas Downsample - sd4 filtrada')
% subplot(2,2,4)
% stem(fd,abs(fft(sdf4)))
% legend('Filter4 (HighPass)')

%% Compand
Mu = 255; % Parameter for mu-law compander

% Segment1
V = max(sdf1);
% 1. Quantize using equal-length intervals and no compander.
[index,quants1,distor] = quantiz(sdf1,0:floor(V),0:ceil(V));

% 2. Use same partition and codebook, but compress
% before quantizing and expand afterwards.
compsig = compand(sdf1,Mu,V,'mu/compressor');
[index,quants1] = quantiz(compsig,0:floor(V),0:ceil(V));
newsig = compand(quants1,Mu,max(quants1),'mu/expander');
newsig = transpose(newsig);
distor2 = sum((newsig-sdf1).^2)/length(sdf1);
% Display both mean square distortions.
[distor, distor2] 

figure();
plot(compsig,':'); % Plot companded signal.
hold on;
plot(sdf1); % Plot original signal.
legend('Companded','Original','Location','NorthWest')

% Segment2
V = max(sdf2);
% 1. Quantize using equal-length intervals and no compander.
[index,quants2,distor] = quantiz(sdf2,0:floor(V),0:ceil(V));

% 2. Use same partition and codebook, but compress
% before quantizing and expand afterwards.
compsig = compand(sdf2,Mu,V,'mu/compressor');
[index,quants2] = quantiz(compsig,0:floor(V),0:ceil(V));
newsig = compand(quants2,Mu,max(quants2),'mu/expander');
newsig = transpose(newsig);
distor2 = sum((newsig-sdf2).^2)/length(sdf2);
% Display both mean square distortions.
[distor, distor2] 

% figure();
% plot(compsig,':'); % Plot companded signal.
% hold on;
% plot(sdf2); % Plot original signal.
% legend('Companded','Original','Location','NorthWest')

% Segment3
V = max(sdf3);
% 1. Quantize using equal-length intervals and no compander.
[index,quants3,distor] = quantiz(sdf3,0:floor(V),0:ceil(V));

% 2. Use same partition and codebook, but compress
% before quantizing and expand afterwards.
compsig = compand(sdf3,Mu,V,'mu/compressor');
[index,quants3] = quantiz(compsig,0:floor(V),0:ceil(V));
newsig = compand(quants3,Mu,max(quants3),'mu/expander');
newsig = transpose(newsig);
distor2 = sum((newsig-sdf3).^2)/length(sdf3);
% Display both mean square distortions.
[distor, distor2] 

% figure();
% plot(compsig,':'); % Plot companded signal.
% hold on;
% plot(sdf3); % Plot original signal.
% legend('Companded','Original','Location','NorthWest')

% Segment4
V = max(sdf4);
% 1. Quantize using equal-length intervals and no compander.
[index,quants4,distor] = quantiz(sdf4,0:floor(V),0:ceil(V));

% 2. Use same partition and codebook, but compress
% before quantizing and expand afterwards.
compsig = compand(sdf4,Mu,V,'mu/compressor');
[index,quants4] = quantiz(compsig,0:floor(V),0:ceil(V));
newsig = compand(quants4,Mu,max(quants4),'mu/expander');
newsig = transpose(newsig);
distor2 = sum((newsig-sdf4).^2)/length(sdf4);
% Display both mean square distortions.
[distor, distor2] 

% figure();
% plot(compsig,':'); % Plot companded signal.
% hold on;
% plot(sdf4); % Plot original signal.
% legend('Companded','Original','Location','NorthWest')


                                    %Mu = 255;
                                    %c1 = compand(sdf1,Mu,max(sdf1),'mu/compressor');
                                    %c2 = compand(sdf2,Mu,max(sdf2),'mu/compressor');
                                    %c3 = compand(sdf3,Mu,max(sdf3),'mu/compressor');
                                    %c4 = compand(sdf4,Mu,max(sdf4),'mu/compressor');

                                    %partition = 0:2^6-1; 
                                    %codebook = 0:2^6;
                                    %[~,~,distor] = quantiz(sd1,partition,codebook);
                                    %[~,quants1] = quantiz(c1,partition,codebook);
                                    %[~,quants2] = quantiz(c2,partition,codebook);
                                    %[~,quants3] = quantiz(c3,partition,codebook);
                                    %[~,quants4] = quantiz(c4,partition,codebook);

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
fComp = inFileStatus.bytes./outFileStatus.bytes