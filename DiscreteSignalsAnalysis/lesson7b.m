%% Lesson 7b = 2021/08/24; STA - Investigating spike-field coupling

% setup
clear, clf, clc
load('SpkBuz.mat')
load('LFPBuz.mat')

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(LFP)); % in s

% plot Rastergram of all neurons
h1 = subplot(3,1,[1 2]);
    for cell=1:40
    spktimes = Raster{cell};
    plot(spktimes,cell*ones(size(spktimes)),'ko')
    hold on
    end
    hold off
    ylabel('Neuron #')

% Plot LFP that identified the 40 neurons
h2 = subplot(3,1,3);
    plot(t,LFP,'k-')

linkaxes([h1 h2],'x')

%% Identify main freq band

% Visual inspection of main freq band
subplot(111)
    [Pxx, F] = pwelch(LFP,4*srate,[],0:0.1:20,srate);
    plot(F,Pxx,'k-')
    ylabel('Power')
    xlabel('Freq (Hz)')

% band identified: 6 to 10 Hz
% filter in this band
theta = eegfilt(LFP',srate,6,10);
% get insantaneous phase
thetaphase = angle(hilbert(theta));

%% Plot a cell spikes over filtred LFP and over extracted phases

spktimes = Raster{7};
spkind = round(spktimes*srate);

h1 = subplot(211);
    plot(t,LFP,'b'); hold on
    plot(t,theta,'r','linew',2)
    plot(t(spkind),theta(spkind),'ko','markerf','g'); hold off
    ylabel('\muV')
    xlabel('Time (s)')

h2 = subplot(212);
    plot(t,rad2deg(thetaphase)+180,'r.'); hold on
    plot(t(spkind),rad2deg(thetaphase(spkind))+180,'ko','markerf','g'); hold off
    ylabel('\Theta Phase (^o)')
    xlabel('Time (s)')

linkaxes([h1 h2],'x')
set(gca,'ytick',0:90:360)

%% Distribution of spikes per theta cycle

clf
spktimes = Raster{7};
spkind = round(spktimes*srate);
spkphase = rad2deg(thetaphase(spkind))+180;
spkphaseRad = thetaphase(spkind)+pi;

% spike-phase histogram
[spkcounts, phasebins]=hist(spkphase,10:20:360);

subplot(221)
    bar([phasebins phasebins+360],[spkcounts spkcounts],'k')
    set(gca,'xtick',0:90:720)
    xlabel('\Theta Phase (^o)')
    ylabel('Spike counts')
    box off

subplot(223)
    rose(spkphaseRad,18)

subplot(222)
    % unitvector = exp(i*2*pi)
    unitvectors = exp(1i*spkphaseRad);
    % subamostrando para plotar
    Iamostra = randi(length(unitvectors),[1,100]);
    compass(unitvectors(Iamostra))







