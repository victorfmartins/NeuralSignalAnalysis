%% 19/08/2021 - Exam Resolution
%/ ************************************************************************** #
%                                                                             #
%                                              ::::::  ::::::    :::    :::   #
%  Exam 1 Resolution                        :+:     :+:         :+:+:  :+:+:  #
%                                          +:+     +:+         +:+ :+ :+ +:+  #
%  By: vfranco- <victorf.martins@usp.br>  +#+     +#+         +#+   +:+ +#+   #
%                                        +#+     +#+         #+#       #+#    #
%  Created: 2021/08/19 13:11              #+#     #+#       #+#       #+#     #
%  Updated: 2021/08/26 00:26                ######  ###### ###       ###.usp  #
%                                                                             #
% *************************************************************************** #

%% Question 1
% Create a vector with the record timestamp for
% this signal, and plot each channel on a separate
% line. Use proper labels for the axes (eg,
% “Time(s)” vs “Channel #”)

% setup
clear, clc, clf
load('MUA.mat')
srate = 30000; % sampling rate, in Hz (1/s)
dt = 1/srate; % sampling period, in seconds (s)
t = (1:length(MUA))*dt; %time vector in s

% plot every channel one above the other
subplot(111)
    plot(t,MUA(1,:)+0,'k'), hold on
    plot(t,MUA(2,:)+1,'b')
    plot(t,MUA(3,:)+2,'g')
    plot(t,MUA(4,:)+3,'r'), hold off
    ax = gca;
    ax.YTick = [0 1 2 3];
    ax.YTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    title('MUA')
    ylabel('mV')
    xlabel('Time (s)')
    xlim([0 60])

%% Question 2
% Using channel 1, make a routine that identifies
% the triggering times of the multi-unit activity
% using as threshold (threshold) the value of
% 0.18. Plot the channel 1 signal along with tags
% indicating the identified trigger times.

% setup
clear, clc, clf
load('MUA.mat')
srate = 30000; % sampling rate, in Hz (1/s)
dt = 1/srate; % sampling period, in seconds (s)
t = (1:length(MUA))*dt;

threshold = 0.18;

% Findpeaks with built in matlab function
[~,Ispikes1]= findpeaks(MUA(1,:),'MinPeakHeight',threshold);

% plot in read only the points discovered to be a spike
subplot(111)
    plot(t,MUA(1,:)-0,'k'), hold on
    hold on
    plot(t(Ispikes1),MUA(1,Ispikes1),'ko','markerf','r','markersize',10)
    hold off
    axis tight
    title('MUA of Channel One')
    ylabel('mV')
    xlabel('Time (s)')


%% Question 3
% Based on the trigger times identified above,
% create a routine that stores the waveforms of
% the action potentials recorded extracellularly
% for each of the 4 channels (always referring to
% the trigger times identified in channel 1). Use
% a ±0.5 ms window centered on the trigger (ie,
% each trigger will have 31 voltage samples). Then
% create a matrix from the concatenation of the 4
% waveforms (for each trigger time) and plot. Note
% that this matrix will have dimension of Number
% of Identified Spikes X 124 voltage samples.


% setup - From previous Questions
clear, clc, clf
load('MUA.mat')
threshold = 0.18;
[~,Ispikes1]= findpeaks(MUA(1,:),'MinPeakHeight',threshold);

% set the size of variables for claimed speed
waveshape1 = zeros(length(Ispikes1),31);
waveshape2 = zeros(length(Ispikes1),31);
waveshape3 = zeros(length(Ispikes1),31);
waveshape4 = zeros(length(Ispikes1),31);
for j=1:length(Ispikes1)
    % take 15 neighbors to each side of the spine
    % center to compose its shape
    waveshape1(j,:) = MUA(1,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape2(j,:) = MUA(2,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape3(j,:) = MUA(3,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape4(j,:) = MUA(4,Ispikes1(j)-15:Ispikes1(j)+15);
end
waveshape = [waveshape1 waveshape2 waveshape3 waveshape4];

subplot(111)
    plot(waveshape','k')
    ax = gca;
    ax.XTick = [15 45 75 105];
    ax.XTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    title('Spike Waveshapes')
    ylabel('mV')
    xlabel('Waveshapes')
    

%% Question 4
% Using the matrix of waveforms constructed above,
% a) compute its covariance matrix, and b) the
% principal components and associated variances,
% ordering them from largest to smallest variance.
% Also, c) plot the first two components (eg,
% using a stem plot), and the scree plot
% (component number vs associated variance) for
% the first 10 principal components. 

% setup - From previous Questions
clear, clc, clf
load('MUA.mat')
threshold = 0.18;
[~,Ispikes1]= findpeaks(MUA(1,:),'MinPeakHeight',threshold);
for j=1:length(Ispikes1)
    waveshape1(j,:) = MUA(1,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape2(j,:) = MUA(2,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape3(j,:) = MUA(3,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape4(j,:) = MUA(4,Ispikes1(j)-15:Ispikes1(j)+15);
end
waveshape = [waveshape1 waveshape2 waveshape3 waveshape4];

% take the covariance matrix of the concatenated
% waveshapes and extract its eigenvectors and
% eigenvalues
C = cov(waveshape);
[PC,lambda] = eig(C);

% Sort from lower to higher variance values to match PC indexes
[variancias, I]=sort(diag(lambda),'descend');
PC1 = PC(:,I(1));
PC2 = PC(:,I(2));

% Seeing PC1 and PC2
clf
subplot(411)
    plot(waveshape','k')
    ax = gca;
    ax.XTick = [15 45 75 105];
    ax.XTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    title('Spike Waveshapes')
    ylabel('mV')
subplot(412)
    stem(PC1)
    ax = gca;
    ax.XTick = [15 45 75 105];
    ax.XTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    title('PC1')
subplot(413)
    stem(PC2)
    ax = gca;
    ax.XTick = [15 45 75 105];
    ax.XTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    title('PC2')
subplot(414)
    scree = variancias(1:10)/sum(variancias);
    plot(100*scree,'b-o','markerfacecolor','b')
    xlabel('Component Number')
    ylabel('% Variance')
    set(gca,'xtick',1:10)
    xlim([0 10+1])
    title('Scree')
% realize that PC1 alone account for almost
% 80% of the varience in the waveshapes
    
%% Question 5
% Then project and plot on a 2D graph the wave
% attributes using the first two principal
% components. Also, use a k-means algorithm to
% identify 2 clusters in this graph.

% setup - From previous Questions
clear, clc, clf
load('MUA.mat')
threshold = 0.18;
[~,Ispikes1]= findpeaks(MUA(1,:),'MinPeakHeight',threshold);
for j=1:length(Ispikes1)
    waveshape1(j,:) = MUA(1,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape2(j,:) = MUA(2,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape3(j,:) = MUA(3,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape4(j,:) = MUA(4,Ispikes1(j)-15:Ispikes1(j)+15);
end
waveshape = [waveshape1 waveshape2 waveshape3 waveshape4];
C = cov(waveshape);
[PC,lambda] = eig(C);
[~, I]=sort(diag(lambda),'descend');
PC1 = PC(:,I(1));
PC2 = PC(:,I(2));


% gather the most and second-most able features to
% differentiate the waveshapes - from the
% perspective that gives each set of points its
% maximum variance
Feature1 = waveshape*PC1;
Feature2 = waveshape*PC2;

numk = 2; % number of clusters
[idx, K] = kmeans([Feature1,Feature2],numk); % cluster algorithm

% Seeing the clusters
subplot(311)
    plot(waveshape','k')
    ax = gca;
    ax.XTick = [15 45 75 105];
    ax.XTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    title('Spike Waveshapes')
    ylabel('mV')
subplot(3,1,[2 3])
    plot(Feature1,Feature2,'ko'); hold on
    plot(Feature1(idx==1),Feature2(idx==1),'ro')
    plot(Feature1(idx==2),Feature2(idx==2),'bo')
    plot(K(1,1),K(1,2),'kx','markersize',20,'linew', 3)
    plot(K(2,1),K(2,2),'kx','markersize',20,'linew', 3); hold off
    xlabel('Feature1')
    ylabel('Feature2')
    title('PC1 vs PC2 Scatter-Plot')
    
    
%% Question 6
% Finally, in one subplot (2,1,1) plot the average
% waveform for each of the clusters (ie, for each
% of the two identified neurons), and in another
% subplot (2,1,2), plot the rastergram of spikes
% for these two cells. For the last graph, use
% “Time (s)” and “Neuron #” as labels. 

% setup - From previous Questions
clear, clc, clf
load('MUA.mat')
srate = 30000;
dt = 1/srate;
t = (1:length(MUA))*dt;
threshold = 0.18;
[~,Ispikes1]= findpeaks(MUA(1,:),'MinPeakHeight',threshold);
for j=1:length(Ispikes1)
    waveshape1(j,:) = MUA(1,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape2(j,:) = MUA(2,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape3(j,:) = MUA(3,Ispikes1(j)-15:Ispikes1(j)+15);
    waveshape4(j,:) = MUA(4,Ispikes1(j)-15:Ispikes1(j)+15);
end
waveshape = [waveshape1 waveshape2 waveshape3 waveshape4];
C = cov(waveshape);
[PC,lambda] = eig(C);
[variancias, I]=sort(diag(lambda),'descend');
PC1 = PC(:,I(1));
PC2 = PC(:,I(2));
Feature1 = waveshape*PC1;
Feature2 = waveshape*PC2;
numk = 2;
[idx, K] = kmeans([Feature1,Feature2],numk);

% colorly separate the different neuron wave
% shapes and give each set a mean
subplot(211)
    plot(waveshape(idx==1,:)','r-'); hold on
    plot(waveshape(idx==2,:)','b-')
    plot(mean(waveshape(idx==1,:))','k-','linew',3)
    plot(mean(waveshape(idx==2,:))','k-','linew',3); hold off
    ax = gca;
    ax.XTick = [15 45 75 105];
    ax.XTickLabel = {'Channel 1','Channel 2','Channel 3','Channel 4'};
    axis tight
    ylabel('mV')
    title('Spike Waveshapes')

subplot(212)
    plot(t(Ispikes1(idx==1)),1*ones(size(Ispikes1((idx==1)))),'ro','markerfacecolor','w')
    hold on
    plot(t(Ispikes1(idx==2)),2*ones(size(Ispikes1((idx==2)))),'bo','markerfacecolor','w')
    hold off
    ax = gca;
    ax.YTick = [1 2];
    ax.YTickLabel = {'Neuron 1','Neuron 2'};    
    xlabel('Time (s)')
    ylim([0.5 2.5])
    title('Rastergram')



    
%% Question 7 - Data set 2
% Make code that plots a rastergram of spikes
% for a single cell, a single taste, across all
% trials.  

% setup - data we want
clear, clc, clf
load('SpikeGustatoryCortex.mat')
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;
taste = 1;

% gather and plot the spike train of each trial
% one above the other
subplot(111)
    count = 0;
    for trial = 1:10
        count = count+1;
        spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 
        plot(spks,(trial+(taste-1)*10)*ones(size(spks)),...
        'ko','markerfacecolor', 'k')
        hold on
    end
    plot([0 0],[0 count+1],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    ylim([0 count])
    xlabel('Time (s)')
    ylabel('Trial #')
    set(gca,'fontsize',14)
    title('Rastergram')


%% Question 8
% Expand the routine to plot cell response to suit
% all tastes. Use some visual aid to indicate
% different tastes (for example, the color of the
% shot mark, or lines separating the trials of
% different tastes).

clear, clc, clf
load('SpikeGustatoryCortex.mat')

% setup - data we want
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;

% set marker of taste
cor{1} = 'c';
cor{2} = 'g';
cor{3} = 'r';
cor{4} = 'b';

% gather and plot the spike train of each trial
% one above the other
subplot(111)
    count = 0;
    for taste = 1:4
    for trial = 1:10
        count = count+1;
        spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 
        plot(spks,(trial+(taste-1)*10)*ones(size(spks)),...
        'ko','markerfacecolor',cor{taste})
        hold on
    end
    end
    plot([0 0],[0 count+1],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    ylim([0 count])
    xlabel('Time (s)')
    ylabel('Trial #')
    set(gca,'fontsize',14)
    title('Rastergram')


%% Question 9
% Create a routine that computes and plots a
% discrete PSTH with the spike counts. Adapt the
% previous routine to plot in a subplot(3,1,[1 2])
% the rastergram, and in the subplot(3,1,3) the
% PSTH.

clear, clc, clf
load('SpikeGustatoryCortex.mat')

% setup - data we want
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;

% set marker of taste
cor{1} = 'c';
cor{2} = 'g';
cor{3} = 'r';
cor{4} = 'b';

% store all spk times by concatenation of vectors
SPKSALL = [];

subplot(3,1,[1 2])
    count = 0;
    for taste = 1:4
    for trial = 1:10
        count = count+1;
        spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 
        plot(spks,(trial+(taste-1)*10)*ones(size(spks)),...
        'ko','markerfacecolor',cor{taste})
        hold on
        % cancatanate spk times
        SPKSALL = [SPKSALL,spks]; 
    end
    end
    plot([0 0],[0 count+1],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    ylim([0 count])
    ylabel('Trial #')
    title('Rastergram')

    
%%% Discrete PSTH
winlength = 0.1;
clear PSTH T
nwintotal = 3.5/winlength;
% Loop over windows and count the number of spikes
% in each
for nwin = 1:nwintotal
    LimE = -1 + (nwin-1)*winlength;
    LimD = -1+winlength + (nwin-1)*winlength; 
    PSTH(nwin) = sum(SPKSALL>LimE & SPKSALL<=LimD);
    T(nwin) = (LimE+LimD)/2;
end

subplot(3,1,3)
    % transforming count in Firing Rate (in Hz)
    bar(T,PSTH/(winlength*count),'k')
    xlabel('Time (s)')
    ylabel('Average Firing Rate (Hz)')
    title('PSTH')
    
    
%% Question 10
% Repeat the above procedure, but using a
% continuous PSTH.

clear, clc, clf
load('SpikeGustatoryCortex.mat')

% setup - data we want
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;

% set marker of taste
cor{1} = 'c';
cor{2} = 'g';
cor{3} = 'r';
cor{4} = 'b';

% store all spk times by concatanation of vectors
SPKSALL = [];

subplot(3,1,[1 2])
    count = 0;
    for taste = 1:4
    for trial = 1:10
        count = count+1;
        spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 
        plot(spks,(trial+(taste-1)*10)*ones(size(spks)),...
        'ko','markerfacecolor',cor{taste})
        hold on
        % cancatanate spk times
        SPKSALL = [SPKSALL,spks]; 
    end
    end
    plot([0 0],[0 count+1],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    ylim([0 count])
    ylabel('Trial #')
    title('Rastergram')

clear PSTH T
winlength = 0.1;
%%% This makes the continuous aspect. it has to be step<=winlength. 
%%% If step=winlength then this plot is discrete rather then continuos
%%% If step<windlength then the continuous aspect starts to be seem
%%% try step = 0.01 to see continuous aspect. 
%%% try step<0.01 to see that nothing happens
step = 0.01;
Nwin = (3.5-winlength)/step+1;

for nwin = 1:Nwin
   LimE = -1 + (nwin-1)*step;
   LimD = LimE + winlength;
   PSTH(nwin) = sum(SPKSALL>LimE & SPKSALL<=LimD);
   T(nwin) = (LimE+LimD)/2;   
end

% transforming count in Firing Rate (in Hz)
FR = PSTH/(winlength*count);
subplot(3,1,3)
    % plot(T,FR)
    fill([-1,T,2.5],[0,FR,0],'c')
    xlabel('Time (s)')
    ylabel('Average Firing Rate (Hz)')
    title('PSTH')
    
    
%% Question 11
% Plot a continuous PSTH using the entire dataset;
% that is, use all firing times of all neurons
% found in the different rats/sessions/hemispheres
% taking into account all trials and all tastes
% (Na,Qu,Su,Ca). Interpret the result and estimate
% the population response latency.

clear, clc, clf
load('SpikeGustatoryCortex.mat')

winlength = 0.1;
step = 0.001;
Nwin = (3.5-winlength)/step+1;

% get all spikes without loop
SPKSALL = [spiketimes{:,:,:,:,:,:}]-1;

for nwin = 1:Nwin
   LimE = -1 + (nwin-1)*step;
   LimD = LimE + winlength;
   PSTH(nwin) = sum(SPKSALL>LimE & SPKSALL<=LimD);
   T(nwin) = (LimE+LimD)/2;
end

% count = hemi*rat*day*taste*neuron*trial
count = 2*3*4*4*8*10;
FR = PSTH/(winlength*count);
subplot(1,1,1)
    fill([-1,T,2.5],[0,FR,0],'c')
    set(gca,'fontsize',14)
    xlabel('Time (s)')
    ylabel('Average Firing Rate (Hz)')
    title('PSTH')
    
    
%% Question 12
% For individual neurons, calculate the mean
% firing rate during the period before taste
% stimuli administration (-1 to 0 seconds),
% immediately after taste administration (0 to +1
% seconds), and later after taste stimuli (+1.5 to
% +2.5 seconds). Express these three spike
% frequencies using a 3-bar graph, one for each
% analyzed period.

clear, clc, clf
load('SpikeGustatoryCortex.mat')

% setup - data we want
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;

% get all spikes count the amouth that fall in
% each time bin [-1 0], [0 1], and [15. 2.5]
SPKSALL = [spiketimes{hemi,rat,session,:,unit,:}]-1;
PSTHpre = sum(SPKSALL>-1  & SPKSALL<=0);
PSTHimm = sum(SPKSALL>0   & SPKSALL<=1);
PSTHpos = sum(SPKSALL>1.5 & SPKSALL<=2.5);

subplot(1,1,1)
    bar(-0.5,PSTHpre/(40),  'k'); hold on
    bar( 0.5,PSTHimm/(40),  'k')
    bar(2,   PSTHpos/(40),'k'); hold off
    xlim([-1.1 2.6])
    set(gca,'fontsize',14)
    xlabel('Time (s)')
    ylabel('Average Spike Counts')
    title('PSTH')


%% Question 13
% From the above results, define the neurons as
% being of type PRE, POST1 and POST2. That is, PRE
% neurons are those that fired the most during the
% 1 second interval preceding the administration
% of the tastes, POST1 neurons fired more than 0
% to +1 s, and POST2 neurons fired more than +1.5
% to 2.5 s. Build and plot 3 PSTHs, one for each
% type of neuron; each PSTH must take into account
% all neurons classified as being of the same
% type.

clear, clc, clf
load('SpikeGustatoryCortex.mat')

% setup - data we want
%
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;

cont_pre = 0;
cont_imm = 0;
cont_pos = 0;

PSTHpreAll = [];
PSTHimmAll = [];
PSTHposAll = [];

for hemi = 1:2
for rat = 1:3
for sess = 1:4
for unit = 1:8
    % count the number of spikes that fell in each time bin
    spks = [spiketimes{hemi,rat,session,:,unit,:}]-1;
    spks_count_pre = sum(spks>-1   & spks<=0);
    spks_count_imm = sum(spks>0    & spks<=1);
    spks_count_pos = sum(spks>1.5  & spks<=2.5);
    
    % get the index order of the bin with the most spikes
    [m, I] = max([spks_count_pre spks_count_imm spks_count_pos]);
    
    %%%% ATENTION: the next if excludes all cells that %%%%
    %%%% does not fire even once!                      %%%%
    if (m~=0)
        % uses the last index to classify the neuron spike times
        if (I==1)
            cont_pre = cont_pre + 1;
            PSTHpreAll = [PSTHpreAll, spks];
        end
        if(I==2)
            cont_imm = cont_imm + 1;
            PSTHimmAll = [PSTHimmAll, spks];
        end
        if(I==3)
            cont_pos = cont_pos + 1;
            PSTHposAll = [PSTHposAll, spks];
        end
    end
end
end
end
end

% make histogram-like counting for each cell type
winlength = 0.1;
bincenters = -0.95:winlength:2.45;
PSTHpre = hist(PSTHpreAll,bincenters);
PSTHimm = hist(PSTHimmAll,bincenters);
PSTHpos = hist(PSTHposAll,bincenters);

subplot(3,1,1)
    bar(bincenters,(PSTHpre)/(1*cont_pre),  'k');
    ylim([0 50])
    xlim([-1 2.5])
    title(['PSTH of ' num2str(cont_pre) ' PRE cells'])
subplot(3,1,2)
    bar(bincenters,(PSTHimm)/(1*cont_imm),  'k')
    ylim([0 50])
    xlim([-1 2.5])
    ylabel('Average Spike Counts')
    title(['PSTH of ' num2str(cont_imm) ' POST1 cells'])
subplot(3,1,3)
    bar(bincenters,(PSTHpos)/(1*cont_pos),'k');
    ylim([0 50])
    xlim([-1 2.5])
    xlabel('Time (s)')
    title(['PSTH of ' num2str(cont_pos) ' POST2 cells'])


%% Question 14 - Data set 3
% Plot a bar graph showing the number of spikes for
% each cell (i.e., “Neuron #” vs “Spike Count”).

clear, clc, clf
load('SpkBuz.mat')

% concatenate the spike count for each neuron
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end

subplot(1,1,1)
    bar(1:length(spks),spks,'k')
    set(gca,'fontsize',14)
    xlabel('Neuron #')
    ylabel('Spike Counts')
    axis tight
    title('Spikes per Neuron')


%% Question 15
% Compute and plot the Inter-spike interval
% histogram (ISI histogram), using 5 ms bin, for
% the 9 neurons with more than 1000 spikes. Use a
% sub-panel for each neuron (eg, subplot(3,3,X)).
% Use 250 ms as the X axis limit.

% setup - From previous Questions
clear, clc, clf
load('SpkBuz.mat')
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end

% get the indexes of those neurons that fire more
% than 1000 times
I = find(spks>1000);

% for each of those neurons compute the time
% between successive firings and plot its
% histogram
for i = 1:9
subplot(3,3,i)
    spks = Raster{I(i)};
    ISI = diff(spks);
    [Counts, ISIbin] = hist(ISI,0:0.005:1);
    bar(ISIbin*1000,Counts,'k')
    if (i>6)
    xlabel('ISI (ms)')
    end
    if (mod(i,3)-1==0)
    ylabel('Counts')
    end
    xlim([0 250])
    set(gca,'fontsize',9)
    title(['Neuron' num2str(I(i)) ': #spikes = ' num2str(sum(Counts)+1)])
end

%% Question 16
% For the same 9 neurons, compute and plot their
% autocorrelograms also in different sub-panels,
% using 500 ms as maximum lag and 10 ms bin.

% setup - From previous Questions
clear, clc, clf
load('SpkBuz.mat')
spks = [];
for ncell = 1:40
spks = [spks length(Raster{ncell})];
end
I = find(spks>1000);


for i = 1:9
% Defining precision and units
binw = 10;                      % in ms
maxlags = 500;                  % in ms
SpkTimes = Raster{I(i)}*1000;   % in ms

% get spike count separated by bin
[spkcount, ~] = hist(SpkTimes,0:binw:(3000*1000));

% compute auto-correlogram with the counted spikes
% of each bin
[ACG, lags] = xcorr(spkcount,maxlags/binw);
ACG(lags == 0) = 0;

subplot(3,3,i)
    bar(lags*10,ACG)
    if (i>6)
    xlabel('Time (ms)')
    end
    if (mod(i,3)-1==0)
    ylabel('Spike Coincidences')
    end
    title(['Neuron' num2str(I(i)) ': #spikes = ' num2str(length(SpkTimes))])
    xlim([-750 750])
    set(gca,'fontsize',9)
end


%% Question 17
% Compute and plot the cross-correlogram between
% neuron 29 (reference) and neuron 36 using 10 ms
% bin and 4 second maximum lag.

clear, clc, clf
load('SpkBuz.mat')

% Defining precision and units
binw = 10;                      % in ms
maxlags = 4000;                 % in ms
SpkTimes29 = Raster{29}*1000;   % in ms
SpkTimes36 = Raster{36}*1000;   % in ms

% get spike count separated by bin
[spkcount29, ~] = hist(SpkTimes29,0:binw:(3000*1000));
[spkcount36, ~] = hist(SpkTimes36,0:binw:(3000*1000));

% compute auto-correlogram with the counted spikes
% of each bin 
[CCG, lags] = xcorr(spkcount36,spkcount29,maxlags/binw);

subplot(111)
    bar(lags*10,CCG)
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    title('Cross Correlogram')
    

%% Question 18
% Get the ISIs of neuron 29, and make a routine to
% scramble the order of the ISIs in order to
% generate “false” trains of fire (surrogate) for
% this neuron. Tip: Use the cumulative sum
% (cumsum) function.

clear, clc, clf
load('SpkBuz.mat')

% compute ISI and shuffle it
spks = Raster{29};
ISI = diff(spks);
ISIr = ISI(randperm(length(ISI)));
[Counts,  ~] = hist(ISI,0:0.005:1);
[Countsr, ISIbin] = hist(ISIr,0:0.005:1);

% compare shuffled and not shuffled ISI
subplot(2,1,1)
    bar(ISIbin*1000, Counts, 'k');  hold on
    plot(ISIbin*1000,Countsr,'ro'); hold off
    xlabel('ISI (ms)')
    ylabel('Counts')
    xlim([0 250])
    title(['Total spikes = ' num2str(sum(Counts))])

% cumulatively sum the interspike interval to
% recreate the spike train. Compare spike train of
% shuffled and not shuffled ISI
subplot(2,1,2)
    plot(cumsum(ISI),2*ones(size(ISI)),'k.'); hold on
    plot(cumsum(ISIr),ones(size(ISIr)),'r.'); hold off
    ax = gca;
    ax.YTick = [1 2];
    ax.YTickLabel = {'Neuron Surrogate','Neuron Actual'};  
    ylim([.5 2.5])
    xlabel('Time (s)')
    ylabel('Neurons')
    title('Rastergram')
    
    
    
%% Question 19 (and 20)
% Loop to compute 100 false trains of triggers for
% neuron 29, and for each compute the 'surrogated'
% cross correlation with neuron 36.

clear, clc, clf
load('SpkBuz.mat')

binw = 10;                      % in ms
maxlags = 4000;                 % in ms

SpkTimes36 = cumsum(diff(Raster{36}*1000)); %in ms
[spkcount36, timebins] = hist(SpkTimes36,0:binw:(3000*1000));

ISI = diff(Raster{29}*1000);

% repeat the ISI shuffling and spike train
% reconstruction 100 times. For each, compute the
% Cross-Correlogram and store it 
subplot(111)
CCGSurrAll = zeros(100,801);
for i=1:100
    SpkTimes29 = cumsum(ISI(randperm(length(ISI))));
    [spkcount29, ~] = hist(SpkTimes29,0:binw:(3000*1000));
    [CCGSurrAll(i,:), ~] = xcorr(spkcount36,spkcount29,maxlags/binw);
end
    
% % Question 20 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, redo question 17, but this time plot
% together, as a statistical threshold, the mean +
% 2 standard deviations of the surrogated
% cross-correlograms.

% Code from Question 17
SpkTimes29 = Raster{29}*1000;   % in ms
SpkTimes36 = Raster{36}*1000;   % in ms
[spkcount29, ~] = hist(SpkTimes29,0:binw:(3000*1000));
[spkcount36, ~] = hist(SpkTimes36,0:binw:(3000*1000));
[CCG, lags] = xcorr(spkcount36,spkcount29,maxlags/binw);

% plot os question 17 plus mean values computed in
% question 19
subplot(111)
    bar(lags*10,CCG); hold on
    plot(lags*10,mean(CCGSurrAll),'r','linew',2)
    plot(lags*10,mean(CCGSurrAll)+2*std(CCGSurrAll),'r','linew',2)
    hold off    
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    title('Cross Correlogram')





















