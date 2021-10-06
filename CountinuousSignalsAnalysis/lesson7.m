%% Lesson 7 = 23/03/2021; Aliasing and Nyquist frequency
% The effect to be observed occurs when changing
% the variable 'plus' above 380 (or when srate2/2
% becomes less than f+plus). What happens with the
% black Pxx is that when the signal's freq is
% higher than the nyquist freq the signal's power
% spectrum has a bounce - the freq shown by the
% power spectrum becomes nyquist - (f-nyquist)

%%% In this Cell we will change the samplingfactor
%%% to understand your PXX 

% Setup
clear, clf, clc
srate = 10000; f = 120; plus = 500;% in Hz
dt = 1/srate; Tmax = 3; t = 0:dt:Tmax; % in Seconds
sinal = sin(2*pi*f*t+pi/2)+sin(2*pi*(f+plus)*t+pi/2);


subplot(211)
    plot(t,sinal)
    xlim([0 0.05])
    xlabel('Time (s)')
    ylabel('mV')
    title('Fake LFP with subsampling')

win = 1*srate;
nfft = 4*srate;
[Pxx, F] = pwelch(sinal,win,[],nfft,srate);

subplot(212)
    plot(F,Pxx,'linew',3)
    xlim([0 800])
    xlabel('Freq (Hz)'),ylabel('Power')
    title('Power Spectrum with Aliasing')

%%% sub-sampling
% when undersampling you reduce the frequency of
% nyquist and if you have any higher frequency in
% your data it will be hit (mirrored) in the new
% frequency of nyquist
samplingfactor = 10;
sinal2 = sinal(1:samplingfactor:end);
srate2 = srate/samplingfactor;
t2 = t(1:samplingfactor:end);

subplot(211)
    hold on
    plot(t2,sinal2,'ko-')
    hold off

win2 = 1*srate2;
nfft2 = 4*srate2;
[Pxx2, F] = pwelch(sinal2,win2,[],nfft2,srate2);

subplot(212)
    hold on
    plot(F,Pxx2,'k-','linew',3)
    hold off
    xlim([0 1500])
    xlabel('Freq (Hz)'),ylabel('Power')

%% In this Cell we will change the freq f to understand your Pxx

% Setup
clear, clf, clc
srate = 10000; % in Hz
dt = 1/srate; Tmax = 3; t = 0:dt:Tmax; % in Seconds

% changing the signal frequency beyond the nyquist
% frequency the result is alising
f = 700;
sinal = sin(2*pi*f*t+pi/2);

subplot(211)
    plot(t,sinal)
    xlim([0 0.05])
    xlabel('Time (s)')
    ylabel('mV')

win = 1*srate;
nfft = 4*srate;
[Pxx, F] = pwelch(sinal,win,[],nfft,srate);

subplot(212)
    plot(F,Pxx,'linew',3)
    xlim([0 800])
    xlabel('Freq (Hz)'),ylabel('Power')
    title('Power Spectrum with Aliasing')

%%% sub-sampling
% changing the sampling factor such that the new
% nyquist frequency is less than my signal
% frequency the result is alaising.
samplingfactor = 10;
sinal2 = sinal(1:samplingfactor:end);
srate2 = srate/samplingfactor;
t2 = t(1:samplingfactor:end);

subplot(211)
    hold on
    plot(t2,sinal2,'ko-')
    hold off
    title(['Sampling rate = ' num2str(srate2) ' Hz'])

win2 = 1*srate2;
nfft2 = 4*srate2;
[Pxx2, F] = pwelch(sinal2,win2,[],nfft2,srate2);

subplot(212)
    hold on
    plot(F,Pxx2,'k-','linew',3)
    hold off
    xlim([0 1500])
    xlabel('Freq (Hz)'),ylabel('Power')

%% Fixing the aliasing effect

% cake recipe: before sub-showing, filter your
% signal below the Nyquist frequency 

% this makes all the higher frequencies that would
% be hit not show up (infecting lower frequencies)
FNyquist = srate2/2;

filtrado = eegfilt(sinal,srate,0,FNyquist);

subplot(211) 
    hold on
    plot(t,filtrado,'r-')
    hold off
    xlim([0.01 0.06])


%%
[Pxx, F] = pwelch(filtrado,win,[],nfft,srate);

subplot(212)
    hold on
    plot(F,Pxx,'r-','linew',3)
    xlim([0 1500])
    xlabel('Freq (Hz)'),ylabel('Power')

% undersampling the filtered signal
samplingfactor = 10;
sinal2 = filtrado(1:samplingfactor:end);
srate2 = srate/samplingfactor;
t2 = t(1:samplingfactor:end);

subplot(211)
    hold on
    plot(t2,sinal2,'ko--')
    hold off

% %
[Pxx2, F] = pwelch(sinal2,win2,[],nfft2,srate2);

subplot(212)
    hold on
    plot(F,Pxx2,'c-','linew',3)
    hold off
    xlim([0 1500])
    xlabel('Freq (Hz)'),ylabel('Power')
