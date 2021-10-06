%% 13/09/2021 - Exam Resolution
%/ ************************************************************************** #
%                                                                             #
%                                              ::::::  ::::::    :::    :::   #
%  Exam 2 Resolution                        :+:     :+:         :+:+:  :+:+:  #
%                                          +:+     +:+         +:+ :+ :+ +:+  #
%  By: vfranco- <victorf.martins@usp.br>  +#+     +#+         +#+   +:+ +#+   #
%                                        +#+     +#+         #+#       #+#    #
%  Created: 2021/09/13 23:11              #+#     #+#       #+#       #+#     #
%  Updated: 2021/09/16 00:26                ######  ###### ###       ###.usp  #
%                                                                             #
% *************************************************************************** #

%% Question 1
% Identify the 9 neurons with more than 1000
% spikes, and compute and plot the STA
% (spike-triggered average) of the theta signal
% for each of them in different sub-panels (eg,
% subplot (3,3,X)). Use 2-second windows centered
% on shots.

% setup
% setup - From previous Questions
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
theta = theta';

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(theta)); % in s

% get number of spikes of each neuron
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end

% get the indexes of those neurons that fire more
% than 1000 times
INeuron = find(spks>1000);

% for each of those neurons compute the time
% between successive firings and plot its
% histogram
for i = 1:9
    spktimes = Raster{INeuron(i)}; % in s
    I = round(spktimes*srate); % round ms to indexes
    Nspike = length(I);
    
    winl = 2*srate;
    STA = zeros(1,winl+1);
    count = 0;
    for nspike=1:Nspike
        % if the given spike is not so much on the
        % edge of the signal that there isn't an
        % entire window around it
    if (I(nspike) > winl/2 & I(nspike)< length(theta) - winl/2)
        count = count+1;
        % take all points around spike that fit to win
        winidx = I(nspike)- round(winl/2):I(nspike)+round(winl/2);
        % only sum all points cause we are only interested in the mean
        STA = STA + theta(winidx)'; % sum, not concatanation
    end
    end
    STA = STA/count;
    
subplot(3,3,i)
    plot(((1:length(STA))-winl/2 )*dt,STA,'b-','linew',3); hold on
    plot([0 0],[min(STA) max(STA)],'k--'); hold off
    if (i>6)
    xlabel('Time (s)')
    end
    if (mod(i,3)-1==0)
    ylabel('mv')
    end
    axis tight
    set(gca,'fontsize',9)
    title(['STA_{Neuron' num2str(INeuron(i)) '}: #spikes = ' num2str(spks(INeuron(i)))])
end


%% Question 2 
% Plot the STA mean of these 9 neurons (Bonus:
% plot additionally ± 1 standard deviation).

% setup
% setup - From previous Questions
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
theta = theta';

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(theta)); % in s

% get number of spikes of each neuron
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end

% get the indexes of those neurons that fire more
% than 1000 times
INeuron = find(spks>1000);

% for each of those neurons compute the time
% between successive firings and plot its
% histogram
STAAll = zeros(9,2*srate+1);
for i = 1:9 %loop over neurons
    spktimes = Raster{INeuron(i)}; % in s
    I = round(spktimes*srate); % round ms to indexes
    Nspike = length(I);
    
    winl = 2*srate;
    STA = zeros(1,winl+1);
    count = 0;
    for nspike=1:Nspike
        % if the given spike is not so much on the
        % edge of the signal that there isn't an
        % entire window around it
    if I(nspike) > winl/2 & I(nspike)< length(theta) - winl/2
        count = count+1;
        % take all points around spike that fit to win
        winidx = I(nspike)- round(winl/2):I(nspike)+round(winl/2);
        % only sum all points cause we are only interested in the mean
        STA = STA + theta(winidx)'; % sum, not concatanation
    end
    end
    STA = STA/count;
    STAAll(i,:) = STA;
end

subplot(111)
    plot(((1:length(STA))-winl/2 )*dt,STAAll,'color',[1 1 1]/1.2); hold on
    plot(((1:length(STA))-winl/2 )*dt,mean(STAAll),'b-','linew',3)
    plot(((1:length(STA))-winl/2 )*dt,mean(STAAll)+1*std(STAAll),'k--')
    plot(((1:length(STA))-winl/2 )*dt,mean(STAAll)-1*std(STAAll),'k--')
    plot([0 0],[min(STA) max(STA)],'k--'); hold off
    xlabel('Time (s)')
    ylabel('mv')
    axis tight
    set(gca,'fontsize',9)
    title(['Mean STA: #spikes = ' num2str(sum(spks(INeuron(:))))])


%% Question 3 
% Using thetaphase, compute and plot the average
% PLV (phase-locking value) for the same 9 neurons
% as a bar graph (Bonus: additionally plot a
% vertical line indicating ± 1 standard
% deviation).

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
thetaphase = thetaphase';

% get number of spikes of each neuron
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end

% get the indexes of those neurons that fire more
% than 1000 times
INeuron = find(spks>1000);

PLVAll = zeros(1,9);
for i = 1:9 %loop over neurons
    spktimes = Raster{INeuron(i)}; % in s
    spkind = round(spktimes*srate); % round ms to indexes
    spkphaseRad = thetaphase(spkind);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    PLVAll(i) = abs(mean(exp(1i*spkphaseRad)));
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

% PLV histogram
[PLVcounts, phasebins]=hist(PLVAll,0:pi/10:2*pi);
M = mean(PLVAll);
SD = std(PLVAll);

subplot(111)
    bar(5,M,'k'); hold on
    errorbar(5,M,SD,'r-sq','markerfacecolor','r'); hold off
    xlim([0 10])
    set(gca,'xtick',[])
    ylabel('Mean PLV')
    title(['Mean PLV of N_{spikes} = ' num2str(length(spktimes))])
    box off


%% Question 4
% Compute the average theta firing phase (=angle)
% for each of the 9 neurons. Plot these 9 values
% using some circular chart (eg, circular histogram [rose
% function in Matlab] or compass [compass
% function]).

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
thetaphase = thetaphase';

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(theta)); % in s

% get number of spikes of each neuron
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end

% get the indexes of those neurons that fire more
% than 1000 times
INeuron = find(spks>1000);

meanAngle = zeros(1,9);
for i = 1:9 %loop over neurons
    spktimes = Raster{INeuron(i)}; % in s
    spkind = round(spktimes*srate); % round ms to indexes
    spkphaseRad = thetaphase(spkind);
    meanAngle(i) = mean(spkphaseRad);
end

% subplot(111)
%     rose(meanAngle,18)
    
subplot(111)
    compass(exp(1i*meanAngle))
%     meanvector = mean(exp(1i*meanAngle)); hold on
%     compass(meanvector,'k'); hold off
    title('Mean Angles')


%% Question 5
% For neuron 7, plot the probability of firings
% per theta phase (use 20 degree bins). Plot on
% the same graph the fit (fitting) of this
% probability by a Von Mises distribution.

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')

% get the spike times and there indexes
spktimes = Raster{7}; % in s
spkind = round(spktimes*srate);  % round ms to indexes

% count the number of spikes in each bin
phasebins = deg2rad(-170:20:170);
[counts] = hist(thetaphase(spkind),phasebins);

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
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('spk probability')
    title(['Kappa_{Von Mises} = ' num2str(kappa)],'fontsize',14)
    
    
%% Question 6
% For neuron 7, compute the “chance” distribution
% for the values of the Von Mises Kappa parameter.
% To do this, compute 1000 surrogated Kappa
% values. Each surrogate can be obtained, for
% example, by using random theta phases as the
% firing phases. Plot the surrogated values
% through a count histogram. Plot on the same
% graph the real Kappa value for neuron 7.

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')

% get the spike times and there indexes
spktimes = Raster{7}; % in s
spkind = round(spktimes*srate); % round ms to indexes
Nspk = length(spktimes);

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


Nsurr = 1000;
% % % % computing the kappa chance keeping ISI structure
% % % ISI = diff(spktimes);
% % % % clear Kappasurr
% % % Kappasurr = zeros(1,Nsurr);
% % % for nsurr = 1:Nsurr
% % %     ISIsurr = ISI(randperm(length(ISI)));
% % %     spksurr = cumsum(ISIsurr);
% % %     spkindsurr = round(spksurr*srate);
% % %     [counts, phasebins] = hist(thetaphase(spkindsurr),phasebins);
% % %     [mu,PLV,sigma,CI,Kappa]=anglemean(thetaphase(spkindsurr));
% % %     Kappasurr(nsurr) = Kappa;
% % % end

% computing the kappa chance (without ISI structure)
% clear Kappasurr
Kappasurr = zeros(1,Nsurr);
for nsurr = 1:Nsurr
    % get Nspk random spikes anywhere in all recording
    spkindsurr = randi(length(thetaphase),[1 Nspk]);
    % take the filtred to theta LFP phases at those indexes
    % compute and store each kappa
    [mu,PLV,sigma,CI,Kappa]=anglemean(thetaphase(spkindsurr));
    Kappasurr(nsurr) = Kappa; 
end

pvalueSurrogado = 1 - sum(kappa > Kappasurr)/length(Kappasurr);

% What is the chance that the realKappa belongs to
% the distribution we got by chance? It is the %
% of the surr that landed higher than the realKappa
subplot(111)
    hist(Kappasurr,0:0.01:1); hold on
    plot([kappa kappa],[0 250],'r','linew',2)
    hold off
    xlim([0 0.5])
    ylabel('Surrogate \kappa Counts')
    xlabel('\kappa')
    set(gca,'fontsize',14)
    title(['statistical_p = ' num2str(pvalueSurrogado) ...
        ' (precision of ' num2str(1/Nsurr) ')'])
    
    
%% Question 7
% Compute, for neuron 10, the entropy-based
% firing-phase coupling index theta.  

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')

% get the spike times and there indexes
spktimes = Raster{10}; % in s
spkind = round(spktimes*srate); % round ms to indexes
Nspk = length(spktimes);

% count the number of spikes in each bin
phasebins = deg2rad(-170:20:170);
[counts, phasebins] = hist(thetaphase(spkind),phasebins);

% scale dow the count to get a prob. distribution
p = counts/sum(counts);

% compute MI and store it
MI = (log(length(phasebins))+sum(p(p>0).*log(p(p>0)))) ...
     /log(length(phasebins));


%% Question 8 
% Do as in question 6, but using neuron 10 and
% entropy-based coupling index. 

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')

% get the spike times and there indexes
spktimes = Raster{10}; % in s
spkind = round(spktimes*srate); % round ms to indexes
Nspk = length(spktimes);

% count the number of spikes in each bin
phasebins = deg2rad(-170:20:170);
[counts, phasebins] = hist(thetaphase(spkind),phasebins);

% scale dow the count to get a prob. distribution
p = counts/sum(counts);

% compute MI and store it
realMI = (log(length(phasebins))+sum(p(p>0).*log(p(p>0)))) ...
     /log(length(phasebins));

Nsurr = 1000;
% % % % computing the MI chance keeping ISI structure
% % % ISI = diff(spktimes);
% % % % clear MIsurr
% % % MIsurr = zeros(1,Nsurr);
% % % for nsurr = 1:Nsurr
% % %     ISIsurr = ISI(randperm(length(ISI)));
% % %     spksurr = cumsum(ISIsurr);
% % %     spkindsurr = round(spksurr*srate);
% % %     [counts, phasebins] = hist(thetaphase(spkindsurr),phasebins);
% % %     p = counts/sum(counts);
% % %     MIsurr(nsurr) = (log(length(phasebins))+sum(p(p>0).*log(p(p>0))))...
% % %                     /log(length(phasebins));
% % % end

% computing the MI chance (without ISI structure)
% (3x faster)
% clear MIsurr
MIsurr = zeros(1,Nsurr);
for nsurr = 1:Nsurr
    % get Nspk random spikes anywhere in all recording
    spkindsurr = randi(length(thetaphase),[1 Nspk]);
    % take the filtred to theta LFP phases at those indexes
    spkphaseRadsurr = thetaphase(spkindsurr);
    % count the number of spikes in each bin
    [counts, phasebins] = hist(spkphaseRadsurr,phasebins);
    % scale dow the count to get a prob. distribution
    p = counts/sum(counts);
    % compute MI and store it
    MIsurr(nsurr) = (log(length(phasebins))+sum(p(p>0).*log(p(p>0))))...
                    /log(length(phasebins));
end

pvalueSurrogado = 1 - sum(realMI > MIsurr)/length(MIsurr);

% What is the chance that the realMI belongs to
% the distribution we got by chance? It is the %
% of the surr that landed higher than the realMI
subplot(111)
    hist(MIsurr,0:0.0001:0.01); hold on
    plot([realMI realMI],[0 75],'r','linew',2)
    hold off
    xlim([0 0.01])
    ylabel('Surrogate MI Counts')
    xlabel('MI')
    title(['statistical_p = ' num2str(pvalueSurrogado) ...
           ' (precision of ' num2str(1/Nsurr) ')'])


%% Question 9 
% Based on the above results, compute the
% percentage of surrogated values (from the
% entropy-based index) that fell below the actual
% measured value.

SurrBelowRealMI = (1 - pvalueSurrogado) * 100

%% Question 10 (Questão Bônus) 
% Using neuron 30 and the PLV metric, investigate
% whether this neuron is more strongly modulated
% by the past, present, or future of the theta
% wave. To do this, compute multiple PLV values
% using lags from –1 second to +1 second. For
% example, the PLV for the –30 ms lag is obtained
% by computing the PLV using the 30 ms theta
% phases before the trigger times. Plot the result
% as a line graph. Also plot a vertical line at
% lag = 0 (ie the PLV value for the present theta
% phase). Identify the lag where the PLV is
% maximum and interpret the result.

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
thetaphase = thetaphase';

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(theta)); % in s

%create time vector
lags = -1:dt:1;

spktimes = Raster{30}; % in s
spkind = round(spktimes*srate); % round ms to indexes

PLVAll = zeros(size(lags));
count = 0;
for lag = lags %loop over lags
    count = count + 1;
    spkphaseRad = thetaphase(spkind+lag*srate);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    PLVAll(count) = abs(mean(exp(1i*spkphaseRad)));
    %%%%%%%%%%%%%%%%%%%%%%%%%
end


subplot(111)
    plot(lags,PLVAll); hold on
    plot([0 0], [0 0.2]); hold off
    xlim([-1 1])
    xlabel('lags (s)')
    ylabel('PLV')
    title(['N_{spikes} = ' num2str(length(spktimes))])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Due to the PLV peak not being at zero lag we can
% ask ourselves if there is a causal relationship
% between the spike and the LFP. As the spike
% comes before the PLV peak, it could be that the
% neuron spike is part of an assembly that
% generates the LFP. This would take a few
% milliseconds to react to the spike stimuli.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Question 11
% Using neuron 36, compute a time series of firing
% rates using 30-second windows with no overlap.
% That is, compute the number of shots between 0
% and 30s, 30s-60s, 60s-90s, etc (and divide each
% count by 30 to get the rate in Hz). Plot the
% result.  

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
thetaphase = thetaphase';

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(theta)); % in s

spktimes = Raster{36};  % in s

timebins = 0:30:max(spktimes);

% count the number of spikes in each bin
[counts, ~] = hist(spktimes,timebins);

FR = counts/30;

subplot(111)
    plot(timebins+15,FR)
    axis tight
    xlabel('Time (s)')
    ylabel('30s Mean Firing Rate (Hz)')
    title(['PSTH N_{spikes} = ' num2str(length(spktimes))])
    
    
%% Question 12 
% Using a subplot(1,2.1), plot the individual
% firing rate values separated into two groups:
% before and after 20 minutes of recording (1200
% seconds).

% setup
clear, clc, clf
load('SpkBuz.mat')
load('thetaBuz.mat')
thetaphase = thetaphase';

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(theta)); % in s

spktimes = Raster{36};  % in s

timebins = 0:30:max(spktimes);

% count the number of spikes in each bin
[counts, ~] = hist(spktimes,timebins);

FR = counts/30;

% get the current state of the random generators.
% every random number generate after a rng(s)
% calling will output the same random values
s = rng; 
for limiar = 0:20;

TP = sum(FR(42:end)>=limiar); % absolute number of true  positives
FN = sum(FR(42:end)< limiar); % absolute number of false negatives

TPR = TP/(FN+TP); % true positive rate (aka sensitivity)

FP = sum(FR(1:41)>=limiar); % false positive
TN = sum(FR(1:41)< limiar); % true negative

FPR = FP/(FP+TN); % false positive rate (aka 1-specificity)

subplot(121)
    % restart the random generator to state s (of
    % before the loop so that the plots dont
    % change in each iteration  
    rng(1);
    % all data lies in a vertical line and many
    % markers can overlap due to the proximity of
    % points. To make them more visible assign a
    % horizontal displacement to each one using
    % randn function
    plot(2*ones(size(FR(42:end)))+ ...
      0.1*randn(size(FR(42:end))),FR(42:end),'rv'); hold on
    plot(1*ones(size(FR(1:41)))+ ...
      0.1*randn(size(FR(1:41))),FR(1:41), 'kv')
    plot([0 3],[limiar limiar ],'m--','linew',3); hold off
    ylim([-2 10])
    ylabel('30s Mean Firing Rate (Hz)')
    xlabel('1200s mark')
    set(gca,'xtick',1:2,'xticklabel',{'Pré';'Pos'},'fontsize',12)
    xlim([0 3])
    title('Labeled Data')
    
subplot(122)
    plot([0 1],[0 1],'k--'); hold on
    plot(FPR,TPR,'ko','markerface','k','markersize',12);
    xlabel('FPR (1-Specificity)')
    ylabel('TPR (Sensitivity)')
    axis square
    set(gca,'fontsize',12)
    title('ROC-Curve points')

pause(0.1)
end
    hold off

%% Question 13
% Next, compute and plot – in subplot(1,2,2) – the
% ROC curve obtained for the classification of
% logging time to be >20 minutes using these
% firing rates. (That is, this ROC curve will tell
% you whether the firing rate of neuron 36 is good
% or not at classifying the recording time as
% being greater or less than 20 minutes).

% get all points, put them in order and use each
% one as treshhold for fine grain precision
threshold = sort(FR);

clear TPR FPR
counts = 0;
for limiar = threshold
    counts = counts+1;
    
TP = sum(FR(42:end)>=limiar); % absolute number of true  positives
FN = sum(FR(42:end)< limiar); % absolute number of false negatives

TPR(counts) = TP/(FN+TP); % true positive rate (aka sensitivity)

FP = sum(FR(1:41)>=limiar); % false positive
TN = sum(FR(1:41)< limiar); % true negative

FPR(counts) = FP/(FP+TN); % false positive rate (aka 1-specificity)  
end

subplot(122)
    hold on
    plot(FPR,TPR,'linew',3)
    hold off
    
    
%% Question 14
% Compute the area under the curve (AUC) from the
% ROC curve obtained above. Use this value as the
% title of the ROC curve graph. 

% AUC = Area Under the Curve
% function trapz makes a trapesoidal (or midd point) integration
realAUC = trapz(FPR(end:-1:1),TPR(end:-1:1));

subplot(122)
    title(['ROC-AUC = ' num2str(realAUC)])
    
    
%% Question 15 (Bônus) 
% Compute and plot 100 surrogated ROC curves, and
% store the AUCs. Tip: note that each firing rate
% is associated with either the label >20 min or
% the label <20 min. To compute a surrogate ROC
% curve, scramble these labels.

threshold = sort(FR);

clear AUCsurr
for i = 1:100
I = (randn(size(FR))>0);
Ilabel1 = FR(find(I==0));
Ilabel2 = FR(find(I==1));

clear TPR FPR
counts = 0;
for limiar = threshold
counts = counts+1;
    
TP = sum(Ilabel1>=limiar);
FN = sum(Ilabel1< limiar);

TPR(counts) = TP/(FN+TP); % true positive rate (aka sensitivity)

FP = sum(Ilabel2>=limiar); % false positive
TN = sum(Ilabel2< limiar); % true negative

FPR(counts) = FP/(FP+TN); % false positive rate (aka 1-specificity)  
end

% AUC = Area Under the Curve
% function trapz makes a trapesoidal (or midd point) integration
AUCsurr(i) = trapz(FPR(end:-1:1),TPR(end:-1:1));
end

% % Question 16 (Bônus) 
% Plote o histograma de contagem de
% valores de AUC-ROC surrogadas junto com o valor
% real medido. Baseado no resultado, podemos dizer
% que o neurônio 36 possui informação sobre o
% tempo de registro ser maior ou menor do que 20
% minutos.

pvalueSurrogado = 1 - sum(realAUC > AUCsurr)/length(AUCsurr)

% What is the chance that the realMI belongs to
% the distribution we got by chance? It is the %
% of the surr that landed higher than the realMI
subplot(111)
    hist(AUCsurr,0:0.01:1); hold on
    plot([realAUC realAUC],[0 10],'r','linew',2)
    hold off
    xlim([0 1])
    ylabel('Surrogate AUC Counts')
    xlabel('AUC-ROC')
    title(['statistical_p = ' num2str(pvalueSurrogado) ...
           ' (precision of ' num2str(1/100) ')'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yes, we can say that neuron 36 has information
% about the recording time being greater or less
% than 20 minutes. Because the chances that he has
% indicated that he has this information by chance
% is quite low, less than 1%.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







































