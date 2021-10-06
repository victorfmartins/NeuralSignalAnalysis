%% Lesson 10 = 08/04/2021; Hilbert Transform
% Hilbert transform
% Takes the fft from the signal and shifts each
% vector of the fft by 90o then takes the ifft  

% Setup
clear, clc, clf, close
srate = 1000;
dt = 1/srate; t = dt:dt:3;

ruidobranco = randn(size(t));
LFP = eegfilt(ruidobranco,srate,8,10);

% The hilbert function, built-in of Matlab,
% outputs the Analytical Signal Representation
% defined by RA = LFP + iH(LFP), where H(LFP) is
% the signal's Hilbert transform

RA = hilbert(LFP);

sinal_original = real(RA);
HT = imag(RA);
Amp = abs(hilbert(LFP));

subplot(1,1,1)
plot(t,LFP)
hold on
% plot(t,sinal_original)
plot(t,HT,'g-')
plot(t,Amp,'k-','linewidth',2)
hold off


%% Visualization of the analytic signal vector and the analytic signal 

for n = 1:2:3000; % amostra
subplot(311)
plot(t,LFP)
hold on
plot(t,HT,'r-')
plot(t,Amp,'k-','linewidth',2)
plot(t(n),LFP(n),'bo')
plot(t(n),HT(n),'ro')
plot(t(n),Amp(n),'ko')
hold off
xlim([0 3])

subplot(3,1,[2 3])
plot([0 LFP(n)],[0 HT(n)],'k-')
hold on
plot(LFP(n),HT(n),'ko')
title(['Amplitude instant?nea = ' num2str(Amp(n))])
axis square
xlim([-.2 .2])
ylim(xlim())
hold off
pause(0.0001)
end


%% #d Representation of the Signal

% A signal enveloped by another
LFP = sin(2*pi*8*t).*sin(2*pi*0.5*t);
HT = imag(hilbert(LFP));
EnvAmp = abs(hilbert(LFP));

subplot(311)
plot(t,LFP)
hold on
plot(t,EnvAmp,'k-','linew',3)
hold off

% Representation 3D
% red line is only the demarcation of zero across time
subplot(3,1,[2 3])
plot3(t,LFP,HT)
hold on
plot3(t,zeros(size(t)),zeros(size(t)),'r-','linew',1)
hold off

xlabel('time (s)')
ylabel('Real')
zlabel('Imag')

%% hilber(constant_signal) = zeros(size(constante_signal))

X = ones(1,1000);
HT = imag(hilbert(X));
AmplitudeEnvelope = abs(hilbert(X));

%% Amplitude on frequencie ff

clear, clc, clf
srate = 1000; ff = 20;
dt = 1/srate; t = dt:dt:10;

LFP = 2*sin(2*pi*20*t) + sin(2*pi*70*t);

order = 300; % ordem of the filter
kernel = sin(2*pi*ff*t(1:order));
kernel = kernel/sum(kernel.^2);

% filters twice to eliminate phase shift 
Filtrado = conv(LFP,kernel,'same');
Filtrado = conv(Filtrado,kernel(end:-1:1),'same');

HT = imag(hilbert(Filtrado));
Amp = abs(hilbert(Filtrado));

% computing the average amplitude (energy measure)
% (most currently used)
MeanAmp = mean(Amp);

% computing the Root Mean Square (energy measure)
% (most used in the past) 
RMS = sqrt(mean(Filtrado.^2));

figure(1)
plot(t,LFP)
hold on
plot(t,Filtrado,'b-')
plot(t,HT,'r-')
plot(t,Amp,'k-','linew',2)
xlim([1 1.2])
hold off

title(['Freq = ' int2str(ff) ...
    ' Hz;  Amp = ' num2str(MeanAmp) ...
    ';  RMS = ' num2str(RMS)])

%% Amplitude Spectrum
% Same as above but with amplitude per frequency

clear, clc, clf
srate = 1000; freqvector = 1:1:100;
dt = 1/srate; t = dt:dt:10;

LFP = sin(2*pi*20*t) + sin(2*pi*40*t)+ sin(2*pi*70*t);

order = 200; % ordem do filtro
for f = freqvector;
kernel = sin(2*pi*f*t(1:order));
kernel = kernel/sum(kernel.^2);

Filtrado = conv(LFP,kernel,'same');
Filtrado = conv(Filtrado,kernel(end:-1:1),'same');

HT = imag(hilbert(Filtrado));
Amp = abs(hilbert(Filtrado));

%%% equal the pwelch (exacly egual?)
MeanAmp = mean(Amp);
AmpSpectrum(f) = MeanAmp; % keep amp at ff

RMS = sqrt(mean(Filtrado.^2));
RMSSpectrum(f) = RMS; % keep amp at ff

subplot(311)
    plot(t,LFP)
    hold on
    plot(t,Filtrado,'b-')
    plot(t,HT,'r-')
    plot(t,Amp,'k-','linew',2)
    xlim([1 1.2])
    hold off

    title(['Freq = ' int2str(f) ...
        ' Hz;  Amp = ' num2str(MeanAmp) ...
        ';  RMS = ' num2str(RMS)])
    pause(0.01)
end

subplot(3,2,[3 5])
    plot(freqvector,AmpSpectrum)
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (mV)')

subplot(3,2,[4 6])
    plot(freqvector,RMSSpectrum)
    xlabel('Frequency (Hz)')
    ylabel('RMS (mV)')


%% Computing the spectrum of the amplitude via eegfilt

clear, clc, clf
srate = 1000; freqvector = 1:1:100;
dt = 1/srate; t = dt:dt:10;
LFP = sin(2*pi*20*t) + sin(2*pi*40*t)+ sin(2*pi*70*t);

% Filtra pass only bandwith frequencies
bandwidth = 4;
clear AmpSpectrum
for flow = freqvector
fhigh = flow+bandwidth;
filtrado = eegfilt(LFP,srate,flow,fhigh);
AmpSpectrum(flow) = mean(abs(hilbert(filtrado)));
end

% the fat diff on the 3 peaks in eegfilt happens
% because the kernel order that is used at each
% freq is equal to 3 cycles of the lowcutoff freq.
% If the frequency increases, each cycle becomes
% tighter in time consuming fewer points = lower
% order. Lower order means less precision in the
% freq domain = higher fat. However, according to
% Perseval's Theorem the area under the curves
% must be equal (there must be an if...) (because
% if it is not false - just look at the picture) 
subplot(211)
plot(freqvector+bandwidth/2,AmpSpectrum)
xlabel('Frequency (Hz)')
ylabel('Amplitude (mV)')

subplot(212)
% Change n at n*srate to change spread
[Pxx, F] = pwelch(LFP,2*srate,[],freqvector,srate);
plot(F,Pxx)
xlim([0 120])


%% TFD "contínous" - time frequency decomposition
% The "non-continuous" PSD computed one PSD per
% window placed each PSD vertically and several
% next to each other. Here we have a measure of
% amplitude over time, so we iterate over freq and
% put the result on top of each other

clear, clc, clf
srate = 1000; freqvector = 1:1:100;
dt = 1/srate; t = dt:dt:10;
LFP = sin(2*pi*10*t) + sin(2*pi*40*t) + sin(2*pi*70*t);
LFP(1:5000) = 0;

clear TFD
tic
count = 0;
for f = freqvector
count = count + 1;

%%% the 1o LFPfilt fix the order at 300 instead of
%%% at 3 cycles
LFPfiltered = eegfilt(LFP,srate,f,f+4,0,300);
% LFPfiltered = eegfilt(LFP,srate,f,f+4);
AmpEnv = abs(hilbert(LFPfiltered));
TFD(count,:) = AmpEnv;
   
% plot(t,LFP)
% hold on
% plot(t,LFPfiltered,'r-')
% plot(t,AmpEnv,'k-','linew',2)
% hold off
% 
% title(['Frequency (Hz) = ' num2str(f+2) ' Hz'])
% xlim([4.5 5.5])
% 
% pause(0.01)
   
end
toc

subplot(211)
imagesc(t,freqvector+2,TFD)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar

% Usando Fourier
subplot(212)
tic
% One must use ridiculously small win 
% to get the same freq resolution as 
% the manual continuous TFD.
% This makes the power scale smaller. Why?
[S,F,T,TFD2]=spectrogram(LFP,0.2*srate,...
            0.18*srate,freqvector+2,srate);
toc
imagesc(t,freqvector+2,TFD2)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar











