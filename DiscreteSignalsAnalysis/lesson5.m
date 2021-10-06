%% Lesson 5 = 2021/08/10; ISI, CCG, ACG
% ISI- Inter-spike interval
% CCG - cross-correlagram
% ACG - auto-correlogram


%% Rastergram
%  Last class we made hastergrams and PSTH like
%  this one. The difference here is the data set
%  that has only 40 cells and have longer
%  recording times. 

% setup
clear, clc, clf
load('SpkBuz.mat')

clear Nspikes
subplot(3,1,[1 2])
    spksall = [];
    for ncell = 1:40
        spks = Raster{ncell};
        Nspikes(ncell) = length(spks);
        plot(spks,ncell*ones(size(spks)),'k.')
        hold on
        spksall = [spksall;spks];
    end
    hold off
    ylabel('Neuron #')
    title('Rastergram')

subplot(313)
    [PSTH, TimeBins] = hist(spksall,5:10:3000);
    bar(TimeBins,PSTH,'k')
    xlabel('Time (s)')
    ylabel('Spike counts')


%% ISI histogram
% cell 27 makes it easier to see two hills in the
% ISI histogram. one around 5ms and other around
% 70ms. The first is the intraburst ISI and the
% second is the interburst ISI
clf
ncell = 27;
spks = Raster{ncell};
% plot(spks,ncell*ones(size(spks)),'k.')

ISI = diff(spks);

[Counts, ISIbin] = hist(ISI,0:0.005:1);

bar(ISIbin*1000,Counts,'k')
xlabel('ISI (ms)')
ylabel('Counts')
xlim([0 250])

title(['Total spikes = ' num2str(sum(Counts))])


%% CCG - cross correlogram

% any two signals
% X = [0 1 1 0 1];
% Y = [1 0 0 1 1];

% any other two signals
% signals need not be the same length
X = [1 2 5 0 12 2 3 1 4 1 32 2 ];
Y = randi(6,[1 20]);

% cross correrogram of X with Y
% The second argument is the one who moves to
% creat the lags
% [CCG, lags]=xcorr(X,Y);

% The thrid argumet gives the maximum number of
% lags (in sampling units)
[CCG, lags]=xcorr(X,Y,10);

subplot(111)
    bar(lags,CCG)
    axis tight
    xlabel('Lag #')
    ylabel('Spike Coincidences')
    title('Cross-Correlogram')

%% ACG  - auto correrogram

% any signal
% unitary signal was chosen to make clear the
% general aspect of auto correlograms of both
% sides decay from lag zero
X = ones(1,10);

% auto-correlagram
[CCG, lags]=xcorr(X);

subplot(111)
    bar(lags,CCG)
    xlabel('Lag #')
    ylabel('Spike Coincidences')
    title('Auto-Correlogram')

%% Spike analysis - looking for syncrony
% in real data any peak in the plot (including in zero lag) 
% indicates some syncrony. 

% create a blank - with no spikes - neuronal
% activity over 50s with sampling rate of 1000Hz
NeuronioA = zeros(1,5000);
NeuronioB = zeros(1,5000);

% randomly get 100 spike moments between 0 and 50s
I = randi(5000,[1 100]);

%%% syncronous case
% assign 1 for TRUE - spike occurence - in the I
% spike times
NeuronioA(I) = 1;
NeuronioB(I) = 1;

IA = find(NeuronioA==1);
IB = find(NeuronioB==1);

subplot(211)
    plot(IA,ones(size(IA)),'k.')
    hold on
    plot(IB,2*ones(size(IB)),'k.')
    hold off
    ylim([0 3])
    xlabel('Time')
    ylabel('Neuron #')
    title('Rastergram')

subplot(212)
    [CCG, lags] = xcorr(NeuronioB,NeuronioA);
    bar(lags,CCG)
    xlim([-100 100])
    ylabel('Coincidences')
    xlabel('Time lag (ms)')
    title('Cross-Correlogram')

    
%% Spike analysis - assyncrony correlation with noise

NeuronioA = zeros(1,5000);
NeuronioB = zeros(1,5000);
I = randi(5000,[1 100]);

%%% assyncronous - where A excites B
NeuronioA(I) = 1;
% No noise leads to replication of last result 20ms shiffited
% Noise represents a shiffit of the positions I by
% a ramdom integer of max 10
noise = randi(10,size(I));
NeuronioB(I+20+1*noise) = 1;
IA = find(NeuronioA==1);
IB = find(NeuronioB==1);

subplot(211)
    plot(IA,ones(size(IA)),'k.')
    hold on
    plot(IB,2*ones(size(IB)),'k.')
    hold off
    ylim([0 3])
    xlabel('Time')
    ylabel('Neuron #')
    title('Rastergram')

subplot(212)
    [CCG, lags] = xcorr(NeuronioB,NeuronioA);
    bar(lags,CCG)
    xlim([-100 100])
    ylabel('Coincidences')
    xlabel('Time lag (ms)')
    title('Cross-Correlogram')

    
%% Infering inibition ralationship between cell with CCG
% Here you have zero syncrony - xcorr gives zero information
% one can see inibition, excitation, oscilational relations
clf
NeuronioA = zeros(1,5000);
NeuronioB = zeros(1,5000);

% 5000 ramdom spikes
NeuronioA(randi(5000,[1 500]))=1;
NeuronioB(randi(5000,[1 500]))=1;

% plus 1000 spikes in neuronA
I = randi(5000,[1,1000]);
NeuronioA(I) = 1;

% each of the 1000 plus spikes leads to inibition
% of neuronB
NeuronioB(I+4) = 0;
NeuronioB(I+5) = 0;
NeuronioB(I+6) = 0;
NeuronioB(I+7) = 0;
NeuronioB(I+8) = 0;
NeuronioB(I+9) = 0;

IA = find(NeuronioA==1);
IB = find(NeuronioB==1);

subplot(211)
    plot(IA,ones(size(IA)),'k.')
    hold on
    plot(IB,2*ones(size(IB)),'k.')
    hold off
    ylim([0 3])
    xlabel('Time')
    ylabel('Neuron #')
    title('Rastergram')

subplot(212)
    [CCG, lags] = xcorr(NeuronioB,NeuronioA);
    bar(lags,CCG)
    xlim([-100 100])
    ylabel('Coincidences')
    xlabel('Time lag (ms)')
    title('Cross-Correlogram')


%% CCG/ACG can be used to see ritmic/oscilatory ralations

%noise plus structure: 1000 random spikes plus
%regular spikes every 20ms
NeuronioA = zeros(1,1000);
NeuronioA(randi(1000,[1 100]))=1;
NeuronioA(20:20:1000) = 1;

IA = find(NeuronioA==1);

subplot(211)
    plot(IA,ones(size(IA)),'k.')
    xlabel('Time')
    ylabel('Neuron #')
    title('Rastergram')

subplot(212)
    [ACG, lags] = xcorr(NeuronioA,200);
    bar(lags,ACG)
    ylabel('Coincidences')
    xlabel('Time lag (ms)')
    title('Auto-Correlogram')

%% Computing ACGs and CCGs in real data

clear, clc, clf
load('SpkBuz.mat')

spks = Raster{40}; % in seconds

% Defining a precision of 1ms
% Does not apply for diferente precisions
SpkTimes = spks*1000; % transform time from seconds into ms
SpkIdx = round(SpkTimes); % make it an integer
disparos(SpkIdx)=1;

maxlags = 500;
[ACG, lags] = xcorr(disparos,maxlags);
% lags == 0 set as true the index in whitch lags = 0.
% making only this point = 0 will plot the graph
% in a more readable maner.
ACG(lags == 0) = 0;
subplot(111)
    bar(lags,ACG,'k')
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    axis tight
    title('Auto-Correlogram')


%% Changing precision

spks = Raster{36}; % in seconds

% Defining precision of 10ms
SpkTimes = spks*1000; % in ms
[spkcount, timebins] = hist(SpkTimes,0:10:(3000*1000));
maxlags = 500/10;
[ACG, lags] = xcorr(spkcount,maxlags);
ACG(lags == 0) = 0;

h1=subplot(311);
    plot(spks*1000,ones(size(spks)),'k.')
    title('Rastergram')
h2=subplot(312);
    plot(timebins,spkcount,'k-')
    title('Histogram')
    linkaxes([h1 h2],'x')
subplot(313)
    bar(lags*10,ACG,'k')
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    axis tight
    title('Auto-Correlogram')

%% Computing CCG:

clear, clc, clf
load('SpkBuz.mat')

% spksA = Raster{7}; % in seconds
% spksB = Raster{36}; % in seconds

% precision of 1 ms
spksA(round(Raster{7}*1000))=1; 
spksB(round(Raster{36}*1000))=1;

subplot(211)
    plot(spksA+1)
    hold on
    plot(spksB)
    hold off
    xlabel('Time (ms)')
    ax = gca;
    ax.YTick = [0 .5 1 1.5 2];
    ax.YTickLabel = {'','Cell 7 ','','Cell 36',''};
    title('Rastergram')

subplot(212)
[CCG, lags] = xcorr(spksB,spksA,2000);
bar(lags,CCG,'k')
xlabel('Time (ms)')
ylabel('Spike Coincidences')
title('Cross Correrogram')

%% Comparing actual CCG with random CCG - Method 1
% the hill in the CCG plot may indicates that cell
% 7 (or 30) excites the other cell or that
% (knowing that this data come from place cells)
% the cell 7 (or 30) comes first in the path of
% the animal.

clear, clc, clf
load('SpkBuz.mat')

% To teste the chance of this hill to happe n by
% chance we need to define what would be by chance
% The chance plot will be the surrogate one
spksA = Raster{7};  % in seconds
spksB = Raster{30}; % in seconds

[PSTH_A, ~] = hist(spksA,0:0.010:3000);
[PSTH_B, timebins] = hist(spksB,0:0.010:3000);

dt = timebins(2)-timebins(1);
maxlag = round(4/dt);
[CCG, lags] = xcorr(PSTH_B,PSTH_A,maxlag);

subplot(211)
    bar(lags*dt*1000,CCG,'k')
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    title('Actual CCG')
    
% Method 1, sorting a number of spikes equal to
% that of the neuron.
% The rigth thing to do hear would be to compute
% this random spikes many times, take the mean and
% the deviation from the mean. For now, lets just
% plot one run of the random algoritm.
% Problem with this approch:
% you destroied all tamporal relations of this
% neuron, thus this is not a good representation
% of any activity of this neuron.
spksB = randi(3000*1000,size(spksB)); % in ms
spksB = spksB/1000; % in seconds
[PSTH_B, ~] = hist(spksB,0:0.010:3000);

subplot(212)
    [CCG, lags] = xcorr(PSTH_B,PSTH_A,maxlag);
    bar(lags*dt*1000,CCG,'k')
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    title('Surrogate CCG')
    ylim([0 60])
    
%% Comparing actual CCG with random CCG - Method 2
% same cells as above

clear, clc, clf
load('SpkBuz.mat')

% To teste the chance of this hill to happen by
% chance we need to define what would be by chance
% The chance plot will be the surrogate one
spksA = Raster{7};  % in seconds
spksB = Raster{30}; % in seconds

[PSTH_A, ~] = hist(spksA,0:0.010:3000);
[PSTH_B, timebins] = hist(spksB,0:0.010:3000);

dt = timebins(2)-timebins(1);
maxlag = round(4/dt);
[CCG, lags] = xcorr(PSTH_B,PSTH_A,maxlag);

subplot(211)
    bar(lags*dt*1000,CCG,'k')
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    title('Actual CCG')
    
% Method 2, circular shift
% does not work in matlab 2016a
[PSTH_B, timebins] = hist(spksB,0:0.010:3000);
PSTH_B = circshift(PSTH_B,100+randi(100*20));

subplot(212)
    [CCG, lags] = xcorr(PSTH_B,PSTH_A,maxlag);
    bar(lags*dt*1000,CCG,'k')
    xlabel('Time (ms)')
    ylabel('Spike Coincidences')
    title('Surrogate CCG')
    ylim([0 60])


%% creating statistical treshold for the CCG

clear CCGSurrAll
for nsurr = 1:100
    nsurr
% Method 1
spksB = randi([3000*1000],[size(spksB)]); % in ms
spksB = spksB/1000; % in seconds
[PSTH_B, timebins] = hist(spksB,0:0.010:3000);

% Method 2
% Only increment position. a good practice would
% be to make it go both ways
% [PSTH_B timebins] = hist(spksB,0:0.010:3000);
% PSTH_B = circshift(PSTH_B,100+randi(100*20));

[CCG, lags] = xcorr(PSTH_B,PSTH_A,maxlag);
CCGSurrAll(nsurr,:) = CCG;
end

subplot(211)
    hold on
    plot(lags*dt*1000,mean(CCGSurrAll),'r','linew',2)
    plot(lags*dt*1000,mean(CCGSurrAll)+2*std(CCGSurrAll),'r','linew',2)
    hold off
    ylim([0 60])






