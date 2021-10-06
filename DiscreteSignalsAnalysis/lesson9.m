%% Lesson 9 = 2021/09/02; Spike-field coupling metrics (continuation)

% set up
clear, clf, clc
load('SpkBuz.mat')
load('LFPBuz.mat')
dt = 1/srate;
t = dt*(1:length(LFP));
theta = eegfilt(LFP',srate,6,10);
thetaphase = angle(hilbert(theta));


%% Von Mises Distribution
% aka the circle's Gaussian (a Gaussiana do c?rculo)

phi = -pi:pi/1000:pi;
% distribution parameters
kappa = 1; % is the inverse of the spread of the distribution
mu = pi/2; % is the mean value of the distribution

% distribution expression
% where besseli is a function that returns a constant
VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));

subplot(111)
    plot(phi,VonMises)
    xlabel('\Phi')
    ylabel('probability')
    xlim([-pi pi])

% check that the sum is 1 (because it is a probability
sum(VonMises*pi/1000)
trapz(phi,VonMises) % function that gives the integral of the distribution

%% 

% get the spike times and there indexes
spktimes = Raster{35};
spkind = round(spktimes*srate);

% count the number of spikes in each bin
phasebins = deg2rad(-170:20:170);
[counts, phasebins] = hist(thetaphase(spkind),phasebins);

% scale dow the count to get a prob. distribution
spkprobability = counts/sum(counts);

% get kappa with tool kit 
[mu,PLV,sigma,CI,kappa]=anglemean(thetaphase(spkind));

% compute the VonMises distribution (yellow in plot)
phi = -pi:pi/1000:pi;
VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));

subplot(111)
    % bar(phasebins,counts,'k')
    bar([phasebins+pi phasebins+pi+2*pi],...
        [spkprobability spkprobability],'k')
    hold on
    plot([phi+pi phi+pi+2*pi],[VonMises VonMises]*deg2rad(20),'y-','linew',3)
    hold off
    xlabel('Phase (rad)')
    ylabel('spk probability')
    title(['Kappa = ' num2str(kappa)],'fontsize',14)


%% index of spike-phase coupling based in entropy
% this metric does not suffer from vector
% aniquilation. The order of the bin count
% distribution in the phases also does not matter
% i.e. if you shuffle the distribution the MI
% continues the same.

% setup (same as above)
spktimes = Raster{35};
spkind = round(spktimes*srate);
phasebins = deg2rad(-170:20:170);
[counts, phasebins] = hist(thetaphase(spkind),phasebins);
p = counts/sum(counts);


% in the entropy information theory 0*log(1) = 0
% so in matlab one has to select only points were
% 0*log(1) does not happen
    % I = p>0
    % H = -sum(p(I).*log(p(I))) % or:
H = -sum(p(p>0).*log(p(p>0)));

% computing maximum entropy
N = length(phasebins);
Hmax = log(N);

% computing the metric
MI = (Hmax-H)/Hmax;

subplot(111)
    bar([phasebins+pi phasebins+pi+2*pi],...
        [p p],'k')
    title(['MI = ' num2str(MI)],'fontsize',14)
    xlabel('Phase (rad)')
    ylabel('spk probability')


%% computando spike-field coupling spectrum

spktimes = Raster{35};
spkind = round(spktimes*srate);

% vector of the frequencies in which PLC and Kappa
% will be taken
freqvector = 0:1:50;

clear Kappaspectrum PLVspectrum MIspectrum
tic
count = 0;
for f = freqvector
    count = count+1

    % filter signal to each freq band and get its phase
    filtrado = eegfilt(LFP',srate,f,f+4);
    phase = angle(hilbert(filtrado));

    %compute PLV and Kappa values
    [mu,PLV,sigma,CI,kappa]=anglemean(phase(spkind));

    % store the PLV and Kappa values
    PLVspectrum(count) = PLV;
    Kappaspectrum(count) = kappa;

    % compute MI and store it
    [counts, phasebins] = hist(phase(spkind),deg2rad(-170:20:170));
    p = counts/sum(counts);
    MI = (log(length(phasebins))+sum(p(p>0).*log(p(p>0))))/log(length(phasebins));
    MIspectrum(count) = MI;
end
toc


%% Plot all 3 metrics for comparissom

subplot(221)
    plot(freqvector+2,PLVspectrum)
    xlabel('Freq (Hz)')
    ylabel('PLV','fontsize',20)

subplot(222)
    plot(freqvector+2,Kappaspectrum)
    xlabel('Freq (Hz)')
    ylabel('\kappa','fontsize',20)

% The proportional distance of maximu point with
% next local minimum is higger in MI than in PLV
% and kappa
subplot(223)
    plot(freqvector+2,MIspectrum)
    xlabel('Freq (Hz)')
    ylabel('MI','fontsize',20)

% One can see that MI goes linear with PLV^2 and kappa^2 
subplot(224)
    plot(Kappaspectrum,MIspectrum,'ko'); hold on
    plot(PLVspectrum,MIspectrum,'ro'); hold off
    xlabel('\kappa (black) and PLV (red)')
    ylabel('MI','fontsize',20)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    