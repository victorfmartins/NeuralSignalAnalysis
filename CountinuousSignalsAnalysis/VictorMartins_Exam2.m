%% 12/0/2021 - Exam Resolution
%/ ************************************************************************** #
%                                                                             #
%                                              ::::::  ::::::    :::    :::   #
%  Exam 2 Resolution                        :+:     :+:         :+:+:  :+:+:  #
%                                          +:+     +:+         +:+ :+ :+ +:+  #
%  By: vfranco- <victorf.martins@usp.br>  +#+     +#+         +#+   +:+ +#+   #
%                                        +#+     +#+         #+#       #+#    #
%  Created: 2021/05/12 13:11              #+#     #+#       #+#       #+#     #
%  Updated: 2021/05/20 00:26                ######  ###### ###       ###.usp  #
%                                                                             #
% *************************************************************************** #

%% Question 1
% (1) (0.05) Calculate and plot the mean evoked
% response (ERP). In the graph, use a vertical
% dashed line to indicate stimulus time.   

clc % clean comm win
clf % clean fig
clear % clean mem
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

ERP = mean(LFP);
m = max(ERP);

plot(t-1, ERP)
hold on
plot([0 0], [-m*1.1 m*1.1], 'k--', 'linew', 2)
ylim([-m*1.1 m*1.1])
hold off

%% Question 2
% (2) (0.15) Compute the average autocorrelogram
% (ACG) (about trials) for the basal period (0-1
% s), for the 1st second after taste
% administration (1-2 s), and for the 2nd second
% (2-3 s). Then plot the average ACG for each of
% these periods in different subplots. Use as
% titles "basal", "early post-taste", and "late
% post-taste"        

clc, clc, clf
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

ACG_basal = zeros(144, 1999);
ACG_early = zeros(144, 1999);
ACG_late = zeros(144, 1999);

for i = 1:144
    disp(i)
    [basal, ~] = xcorr(LFP(i,1:1000));
    ACG_basal(i,:) = basal;

    [early, ~] = xcorr(LFP(i,1001:2000));
    ACG_early(i,:) = early;

    [late, lags] = xcorr(LFP(i,2001:3000));
    ACG_late(i,:) = late;
end

meanACG_basal = mean(ACG_basal);
meanACG_early = mean(ACG_early);
meanACG_late = mean(ACG_late);

subplot(131)
    plot(lags*dt,meanACG_basal)
    xlabel('Lag (ms)')
    ylabel('ACG')
    title('basal')
subplot(132)
    plot(lags*dt,meanACG_early)
    xlabel('Lag (ms)')
    ylabel('ACG')
    title('early post-taste')
subplot(133)
    plot(lags*dt,meanACG_late)
    xlabel('Lag (ms)')
    ylabel('ACG')
    title('late post-taste')

%% Question 3
% (3) (0.1) Compute an ERP time-frequency
% decomposition (TFD) using the continuous wavelet
% transform. Use the Morlet wavelet, and use
% 4:0.2:20 Hz as the frequency vector. Plot the
% result in the subplot (2,1,1). Use a dashed line
% to indicate stimulus timing.     

clc, clc, clf
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

ERP = mean(LFP);

Fc = 1;
NumeroCiclos = 5;
sd = NumeroCiclos/(2*pi*Fc);
Fb = 2*sd^2;
Freq = 4:0.2:20;
scale = Fc./(Freq*dt);
wname = ['cmor' num2str(Fb) '-' num2str(Fc)];

WT = cwt(ERP,scale,wname);
TFD = abs(WT);

imagesc(t-1,Freq,TFD)
hold on
plot([0 0], [0 20], 'k--', 'linew', 2)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(['TFD of ERP with ' num2str(NumeroCiclos) ' cycles wavelet'])
colorbar

%% Question 4
% (4) (0.2) Compute wavelet transforms and TFDs as
% above, but for each trial individually. Then
% plot the average of the TFDs in the subplot
% (2,1,2). Use a dashed line to indicate stimulus
% timing.    

clc, clc, clf
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

ERP = mean(LFP);

Fc = 1;
NumeroCiclos = 5;
sd = NumeroCiclos/(2*pi*Fc);
Fb = 2*sd^2;
Freq = 4:0.2:20;
scale = Fc./(Freq*dt);
wname = ['cmor' num2str(Fb) '-' num2str(Fc)];

TFD = zeros(81, 3500, 144);
for i = 1:144
disp(i)
WT = cwt(LFP(1,:),scale,wname);
TFD(:,:,i) = abs(WT);
end

meanTFD = mean(TFD,3);

WT = cwt(ERP,scale,wname);
ERP_TFD = abs(WT);

subplot (2,1,1)
imagesc(t-1,Freq,ERP_TFD)
hold on
plot([0 0], [0 20], 'k--', 'linew', 2)
hold off
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(['TFD of ERP with ' num2str(NumeroCiclos) ' Cycles Wavelet'])
colorbar

subplot (2,1,2)
imagesc(t-1,Freq,meanTFD)
hold on
plot([0 0], [0 20], 'k--', 'linew', 2)
hold off
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(['Trial Mean TFD with ' num2str(NumeroCiclos) ' Cycles Wavelet'])
colorbar

%% Question 5
% (5) (0.1) Filter each trial between 8 and 12 Hz.
% Plot the average of the filtered signal as well
% as its amplitude envelope. Use a dashed line to
% indicate stimulus timing.    

clc, clc, clf
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

LFP_Filtrado = zeros(144,3500);
for i = 1:144
LFP_Filtrado(i,:) = eegfilt(LFP(i,:),srate,8,12);
end

meanFiltrado = mean(LFP_Filtrado);
Amp = abs(hilbert(meanFiltrado));

m = max(Amp);

plot(t-1,meanFiltrado)
hold on
plot(t-1,Amp)
plot([0 0], [-m*1.1 m*1.1], 'k--', 'linew', 2)
ylim([-m*1.1 m*1.1])
title('Bandpass Filtered (8-12hz) Mean Signal')
hold off

%% Question 6
% (6) (0.1) Filter each trial between 8 and 12 Hz
% and get the amplitude envelope for each trial.
% Plot the average of the amplitude envelopes.  

Amp = zeros(144,3500);
for i = 1:144
Amp(i,:) = abs(hilbert(LFP_Filtrado(i,:)));
end

meanAmp = mean(Amp);
m = max(meanAmp);


%%% The meanAmp amplitude is greater than the
%%% meanFiltered amplitude because in the first
%%% there is no energy cancellation due to phase
%%% misalignment   

% plot(t-1,meanFiltrado)
% hold on
plot(t-1,meanAmp)
hold on
plot([0 0], [-m*1.1 m*1.1], 'k--', 'linew', 2)
ylim([-m*1.1 m*1.1])
title('Bandpass Filtered (8-12hz) Trial Mean Signal')
hold off

%% Bonus Question 1
% (Bonus question 1) (0.1) Plot along the mean
% amplitude ± 1 standard error of the mean (SEM),
% where the SEM is defined as the standard
% deviation (SD) divided by the square root of the
% number of trials.    

clf
SEM = std(Amp)/sqrt(144);

curve1 = meanAmp + SEM;
curve2 = meanAmp - SEM;
x2 = [t-1, fliplr(t-1)];
inBetween = [curve1, fliplr(curve2)];


fill(x2, inBetween, 'cyan');
hold on;
plot(t-1, meanAmp, 'r');
plot([0 0], [-m*1.1 m*1.1], 'k--', 'linew', 2)
ylim([-m*1.1 m*1.1])
title('Bandpass Filtered (8-12hz) Trial Mean Signal with SEM Bounds')
hold off
%% Bonus Question 2
% (Bonus question 2) (0.1) Also compute the
% baseline mean amplitude value as well as its
% standard deviation (SD) using the period from
% 200 ms to 800 ms. Plot the baseline amplitude ±
% 2 SD across continuous horizontal lines.    

clf
SEM = std(Amp)/sqrt(144);

meanBasalAmp = mean(Amp(:,200:800));
Basal = mean(meanBasalAmp);
SDBasal = std(meanBasalAmp);

curve1 = meanAmp + SEM;
curve2 = meanAmp - SEM;
x2 = [t-1, fliplr(t-1)];
inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween, 'cyan');
hold on;
plot(t-1, meanAmp, 'r');
plot([0 0], [-m*1.1 m*1.1], 'k--', 'linew', 2)
plot([0 3.5]-1, [Basal Basal], 'r-', 'linew', 2)
plot([0 3.5]-1, [Basal+2*SDBasal Basal+2*SDBasal], 'y-', 'linew', 2)
plot([0 3.5]-1, [Basal-2*SDBasal Basal-2*SDBasal], 'y-', 'linew', 2)
ylim([-m*1.1 m*1.1])
title('Bandpass Filtered (8-12hz) Trial Mean Signal with SEM Bounds')
hold off

%% Question 7
% (7) (0.15) Filter each trial between 8 and 12
% Hz. Compute and plot in (2,1,1) the
% instantaneous phases of the average of the
% filtered signal, using the first 1.5 seconds
% after the stimulus as the limit of the X axis.
% Then compute the instantaneous frequency of the
% average of the filtered signal and plot in
% (2,1,2), using as the Y-axis limit from 8 to 12
% Hz, and the same X-axis limit as the one in the
% graph above (the first 1.5 seconds after
% stimulus).   

clc, clc, clf
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

LFP_Filtrado = zeros(144,1500);
for i = 1:144
LFP_Filtrado(i,:) = eegfilt(LFP(i,1001:2500),srate,8,12);
end

Phase = angle(hilbert(mean(LFP_Filtrado)));
Freq = angle(exp(1i*diff(Phase)))/(2*pi*dt);

subplot(2,1,1);
plot(t(1001:2500)-1,Phase,'b.')
xlabel('Time (s)')
ylabel('Phase (rad)')

subplot (2,1,2)
plot(t(1001:2499)+dt/2-1,Freq,'b','linew',2)
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% Question 8
% (8) (0.15) Compute and plot the coherence level
% between trials (ITC – inter trial coherence) for
% the 8 to 12 Hz band  


clc, clc, clf
load('GC_LFPs.mat')

srate = 1000;
dt = 1/srate;  % in s
Tmax = length(LFP)/srate;  % in s
t = dt:dt:Tmax;  % in s

LFP_Filtrado = zeros(144,length(t));
PhaseAll = zeros(144,length(t));
for j = 1:144
    LFP_Filtrado(j,:) = eegfilt(LFP(j,:),srate,8,12);
    Phase = angle(hilbert(LFP_Filtrado(j,:)));
    PhaseAll(j,:)=Phase;
end

ITC = abs(mean(exp(1i*PhaseAll)));

plot(t-1,ITC,'k','linew',2)
hold on;
plot([0 0], [0 1], 'r--', 'linew', 2)
xlabel('Time (s)')
ylabel('ITC')
ylim([0 1])

















%%%%%%%%%%%%%%%% End %%%%%%%%%%%%%%%%%%%%%%