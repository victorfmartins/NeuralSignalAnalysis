%% Lesson 8 = 2021/08/31; Spike-field coupling metrics

%%% Class 7b %%%
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
    plot(spktimes,cell*ones(size(spktimes)),'k.')
    hold on
    end
    hold off
    ylabel('Neuron #')
    title('Rastergram')

% Plot LFP that identified the 40 neurons
h2 = subplot(3,1,3);
    plot(t,LFP,'k-')
    title('LFP')

linkaxes([h1 h2],'x')

%% Identify main freq band (class 7b)

% Visual inspection of main freq band
subplot(111)
    [Pxx, F] = pwelch(LFP,4*srate,[],0:0.1:20,srate);
    plot(F,Pxx,'k-')
    ylabel('Power')
    xlabel('Freq (Hz)')
    title('Power Spectrum')

% band identified: 6 to 10 Hz
% filter in this band
theta = eegfilt(LFP',srate,6,10);
% get insantaneous phase
thetaphase = angle(hilbert(theta));

%% Distribution of spikes per theta cycle (class 7b)

clf
spktimes = Raster{1};
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
    title(['N_{spikes} = ' num2str(length(spktimes))])
    box off

subplot(223)
    rose(spkphaseRad,18)

subplot(222)
    % unitvector = exp(i*2*pi)
    unitvectors = exp(1i*spkphaseRad);
    % undersampling to plot (otherwise poluted figure)
    Iamostra = randi(length(unitvectors),[1,100]);
    compass(unitvectors(Iamostra))
    % atention: the mean vector ploted below is
    % always the same no matter the undersample
    meanvector = mean(unitvectors); hold on
    compass(meanvector,'k'); hold off
    title('Subsampled vectors')

%% Phase-Locking Value (PLV)

%%%%%%%%%%%%%%%%%%%%%%%%%
% PLV = abs(meanvector);
PLV = abs(mean(exp(1i*spkphaseRad)));
%%%%%%%%%%%%%%%%%%%%%%%%%

Nspk = length(spkind);
%giving figure proper title
subplot(221)
    title(['PLV = ' num2str(PLV) '; N_{spikes} = ' num2str(Nspk)])
    set(gca,'fontsize',12)

%% Testing PLV Statistical Significance

%%% by Rayleigh test %%%
pvalue = circTestR(thetaphase(spkind));
%giving figure proper title
subplot(223)
    title(['p = ' num2str(pvalue)])

    
%%% by surrogates analysis %%%
% "manual statistical analysis" 
PLV = abs(mean(exp(1i*spkphaseRad)));

clear PLVsurr
% the number of surrogates gives us the porcentage
% precision of our statistical test. e.g. if we
% get p = 0.0074 with Nsurr = 100 then the right p
% is rounded to 0.01. if Nsurr = 1000 then p is
% rounded to 0.007 etc.
Nsurr = 10000;
for nsurr = 1:Nsurr
    % get Nspk random spikes anywhere in all recording
    spkindsurr = randi(length(thetaphase),[1 Nspk]);
    % take the filtred LFP phases at those indexes
    spkphaseRadsurr = thetaphase(spkindsurr)+pi;
    % compute and store each PLV
    PLVsurr(nsurr) = abs(mean(exp(1i*spkphaseRadsurr)));  
end

pvalueSurrogado = 1 - sum(PLV > PLVsurr)/length(PLVsurr)
pvalueRayleigh = circTestR(thetaphase(spkind))

% What is the chance that the realPLV belongs to
% the distribution we got by chance? It is the %
% of the surr that landed higher than the realPLV
subplot(224)
    hist(PLVsurr,0:0.01:1); hold on
    plot([PLV PLV],[0 300],'r','linew',2)
    % annotation('arrow',[PLV PLV],[20 0])
    hold off
    xlim([0 1])
    ylabel('Surrogate PLV Counts')
    xlabel('PLV')
    title(['statistical_p = ' num2str(pvalueSurrogado)])
    

%% Computing the distribution of the preferred spiking phase over neurons

% function that gives many of the values we usualy compute
[PrefSpkPhase,PLV,sigma,confangle,kappa]=anglemean(spkphaseRad);
% e.g. the following two whays give exacly the same
% results than the function above
PLV2 = abs(mean(exp(1i*spkphaseRad)));
PrefSpkPhase2 = angle(mean(exp(1i*spkphaseRad)));

%for all neurons get their significance,
%#ofSpikes, PVL and prefSpkPhase, respectivly.
clear pvalueRayleigh NspkAll SpkPhaseAll
for neuron = 1:40
spktimes = Raster{neuron};
spkind = round(spktimes*srate);
spkphaseRad = thetaphase(spkind)+pi;

pvalueRayleigh(neuron) = circTestR(spkphaseRad);
NspkAll(neuron) = length(spkind);    
PLVAll(neuron) = abs(mean(exp(1i*spkphaseRad)));
SpkPhaseAll(neuron) = angle(mean(exp(1i*spkphaseRad)));
end


%% For every neuron show: pvalue, Nspk, PLV, and PrefSpkPhase

subplot(221)
    bar(NspkAll)
    ylabel('Number of spikes')
    xlabel('Neuron #')
    box off

subplot(222)
    bar(PLVAll)
    ylabel('PLV')
    xlabel('Neuron #')
    box off

subplot(223)
    bar(pvalueRayleigh)
    hold on
    plot([0 41],[0.05 0.05],'k--')
    plot([0 41],[0.01 0.01],'r--')
    ylim([0 0.1])
    hold off
    box off
    ylabel('p-value')
    xlabel('Neuron #')

subplot(224)
    bar(SpkPhaseAll)
    ylabel('Preferred spiking phase')
    xlabel('Neuron #')
    box off

%% Select only significant and relayble neurons for PLV

% selec neurons statisticaly significant for PLV
alpha = 0.01;
I1 = find(pvalueRayleigh<alpha);

% selec neurons with more than 50 spikes
I2 = find(NspkAll>50);

% take intersection of both
I = intersect(I1,I2);

% off those neurons make histogram and compass plot of PrefSpkPhase
clf
subplot(121)
    hist(rad2deg(SpkPhaseAll(I))+180,10:20:350)
    xlabel('Preferred spiking \theta phase (^o)')
    ylabel('Counts')
    box off
    set(gca,'xtick',0:45:360)

subplot(122)
    compass(PLVAll(I).*exp(1i*(SpkPhaseAll(I)+pi)))

% in a deepper level check if the PrefSpkPhase off
% statisticaly significant neurons is itself
% significant
circTestR(SpkPhaseAll(I));

%% PLV Spectrum

spktimes = Raster{35};
spkind = round(spktimes*srate);
freqvector = 0:1:50;

tic
count = 0;
for f = freqvector
count = count+1;
% Since the computation in this loop takes a while
% it was made a test with the following line with
% only 10% of the data (5min) so that one could
% see if the computation was right   
% filtrado = eegfilt(LFP(1:5*60*srate)',srate,f,f+4);
filtrado = eegfilt(LFP',srate,f,f+4);

phase = angle(hilbert(filtrado));
PLV = abs(mean(exp(1i*phase(spkind))));
PLVspectrum(count) = PLV;
end
toc

%

clf
subplot(111)
    % plot(freqvector+2,PLVspectrum,'b-','linew',2)
    % hold on
    plot(freqvector+2,PLVspectrum.^2,'r-','linew',2)
    % hold off
    xlabel('Frequency (Hz)')
    ylabel('Phase-Locking Value')
    ylabel('MRL') % MRL - mean resultant length (other name for PVL)
    ylabel('PLV^2')












