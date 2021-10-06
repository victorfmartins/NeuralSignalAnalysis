%% Lesson 4 = 11/03/2021

%% running LFP

% setup
clear, clf, clc, load('LFP_HG_HFO.mat')
srate = 1000; % in Hz
dt = 1/srate; t = dt:dt:length(lfpHG)*dt; % in s

plot(t,lfpHG)
hold on
plot(t,lfpHFO-1.2)
hold off
box off
xlim([0 1])

% for the animation
stepsize = 0.2;
for nstep = 0:100

xlim([0 1]+nstep*stepsize)
   
ylim([-2 1])
pause(0.1)

end

%% Demonstration of the effect of the non-stationarity of waves
% Here the PSD is not a good tool!
% goal: create a tool that informs us
% on aspects derived from non-seasonality
% of the brain - the TFD

% setup
clear, clf, clc, load('LFP_HG_HFO.mat')
srate = 1000;  f = 8; % in Hz
dt = 1/srate; Tmax = 5; t = dt:dt:Tmax; % in s

% signal
% with amplitude 1 we have a signal present only
% at the end 
LFP = 1*sin(2*pi*f*t); 
LFP(1:2500) = 0;
% with amplitude .7 one signal in every signal
LFP = LFP + 0.7*sin(2*pi*(f+8)*t); 

subplot(311)
    plot(t,LFP)
    xlabel('Time (s)')
    ylabel('mV')
    title('LFP')

% changing the window size influences the analysis
win = 0.5*srate;
overlap = [];
nfft = 10000;

[Pxx, F] = pwelch(LFP,win,overlap,nfft,srate);

% we can see that the power of the two waves is
% equal. This happens in the particular case where
% the shorter existence time of one wave is
% compensated by the smaller amplitude of the
% other. The longer the lifetime, the greater the
% power, the greater the amplitude, the greater
% the power 
subplot(3,1,[2 3])
    plot(F,Pxx,'ko-') 
    xlabel('Frequency (Hz)')
    ylabel('Power')
    xlim([0 30])
    title('Power Spectrum')

%% Time-Frequency Power Decomposition (TFD)
% As in PSD the signal will be split into windows
% and each window will have its sign decomposed by
% intensity of each frequency present.
% The PSD takes the intensities of a frequency
% of all windows and averages to
% give the overall intensity of that frequency.
% TFD does not take this average, it shows in a
% matrix the intensity value in which one axis is
% frequency and the other the windows (time)

% is the factor that determines the trade off freqXtime
window_length = 1*srate;

% is how far a window will go in relation to the
% previous one. USING OVERLAP INCREASES YOUR
% PERCEPTION OF THE TIME WHERE THINGS OCCUR. As
% the analysis is done on the average in relation
% to a window, if a +-localized event happens and
% your window doesn't catch it completely you lose
% this higher power which can be statistically
% significant.
step_size = 0.1*window_length; % (ie., 90% overlap)

% is the number of windows I go
% have (depends on window and step size)
Nwindows = (length(LFP)-window_length)/step_size+1;

clear T TFD
for nwin = 1:Nwindows
    winidx = (1:window_length) + (nwin-1)*step_size;
    [Pxx, F] = pwelch(LFP(winidx),window_length,[],nfft,srate);
    T(nwin) = t(winidx(window_length/2));
    TFD(nwin,:) = Pxx;     
end

subplot(3,1,[2 3])
    imagesc(T,F,TFD')
    axis xy
    ylim([0 30])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colorbar % can uncomment to see color scales
    title('Time-Frequency Power Decomposition (TFD)')
    

%% Show the PSD of each TFD window

for nwin = 1:Nwindows
   plot(F,TFD(nwin,:))
   xlim([0 30])
    title(['Tempo = ' num2str(T(nwin)) ' s'])
   pause(0.1)
end
    

%% using MATLAB's built in function
 
window_length = 1*srate;
% step_size = 0.1*window_length; % (ie., 90% overlap)
overlap = 0.9*window_length;
nfft = 10000;

[~,F,T,TFD]=spectrogram(LFP,window_length,overlap,nfft,srate);

subplot(3,1,[2 3])
    imagesc(T,F,TFD)
    axis xy
    ylim([0 30])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colorbar
    title('Time-Frequency Power Decomposition (TFD)')


%% Normalizations of the Pxx
% when dealing with different animals or different
% situations the absolute value of the parameters
% under study may vary and hence complicates the
% analysis. the normalization makes the parameters
% relative to the general conditions of each
% analysis which makes the comparison between
% analyzes possible.

% be careful with normalization: it may be that
% when offering a certain stimulus the relative
% power of a band changes. This does not mean that
% the power in that band has changed, it may mean
% that the stimulus has changed the power in
% another band which changes the relative power
% between bands.

% this cell shows how changing the power on
% another freq (15Hz) changes the relative power
% of the original freq (8Hz). That's the
% difference between LFP and LFP2
LFP = sin(2*pi*8*t);
LFP2 = sin(2*pi*8*t)+sin(2*pi*15*t);

[Pxx, ~] = pwelch(LFP,2*srate,[],nfft,srate);
[Pxx2, F] = pwelch(LFP2,2*srate,[],nfft,srate);

subplot(2,1,1)
    plot(F,Pxx,'ko-')
    hold on
    plot(F,Pxx2,'ro-')
    hold off
    xlim([0 20])
    ylabel('Power (mV^2/Hz)')
    xlabel('Frequency (Hz)')
    title('Power Spectrum')

% Percentage power
NormPxx = Pxx/sum(Pxx);
NormPxx2 = Pxx2/sum(Pxx2);

subplot(2,1,2)
    plot(F,NormPxx,'ko-')
    hold on
    plot(F,NormPxx2,'ro-')
    hold off
    xlim([0 20])
    ylabel('Normalized Power (%)')
    xlabel('Frequency (Hz)')
    title('Normalized Power Spectrum')


%% Normalizations in TFD by basiline

clear 
srate = 1000;
dt = 1/srate; t=dt:dt:10;

LFP = sin(2*pi*8*t);
LFP(2500:end) =  LFP(2500:end)+0.1*sin(2*pi*50*t(2500:end));

window_length = 1*srate;
overlap = 0.9*window_length;
nfft = 10000;

% TFD has frequencies in the first dimension and
% time in the second
[S,F,T,TFD]=spectrogram(LFP,window_length,overlap,nfft,srate);

% plots with displacement so that the induced
% stimulus is at moment zero 
subplot(211)
    imagesc(T-2.5,F,TFD); hold on
    % marks the moment of stimulation
    plot([0 0],[0 100],'w--','linew',2); hold off
    axis xy
    ylim([0 100])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('TFD')


% baseline normalization
% my event occurred at 2.5s time so i want to know
% the difference before/after 
I = find(T<2.5);

% PSD average over time, that is, every
% contribution to the power of each frequency that
% was made before the stimulus will be deducted
% from the post estimate so that only the
% difference in state is highlighted 
PSDbaseline = mean(TFD(:,I),2);

% takes the average power from before the stimulus
% of every point of PDT 
NormTFD = TFD./repmat(PSDbaseline,1,length(T));

subplot(212)
    imagesc(T-2.5,F,NormTFD); hold on
    plot([0 0],[0 100],'w--','linew',2); hold off
    axis xy
    ylim([0 100])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Baseline Normalized TFD')









