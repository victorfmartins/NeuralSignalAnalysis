%% Lesson 15 = 27/04/2021; Cross and Auto Correlation

% Cross Correlation.
% Can be done with signs of the same size or with
% signs of different size  
clear, clc, clf
X = ones(1,10);
Y = ones(1,3);

% Non-commutative function.
% invert the parameters mirrors the result 
[CCG, lags] = xcorr(X,Y);
% [CCG, lags] = xcorr(Y,X);

% The 3rd argument (optional) of the function says
% the maximum lag to be analyzed  
% [CCG, lags] = xcorr(X,Y,10);

% note that the lag is relative to the time series
% in the second argument of the xcorr() function  

plot(lags,CCG,'ko-')
xlabel('Lag (ms)')
ylabel('CCG')
ylim([0 max(CCG)])

%% Autocorrelograma:

clc, clf
dt = 1/1000;
t = dt:dt:1;

X = ones(1,35);
X = sin(2*pi*5.*t);

% test with 3rd parameter 10 and without it 
% [ACG, lags] = xcorr(X,X,10);
[ACG, lags] = xcorr(X);

% xcorr() accept only one argument 
% [ACG, lags] = xcorr(Y,10);

plot(lags,ACG,'ko-')
xlabel('Lag (ms)')
ylabel('ACG')
ylim([0 max(ACG)])

%% Applications...
% The CCG can be used for:
% 1) check synchrony (when the peak is at zero)
% 2) find causality (when the peak is not at zero)
% 3) check oscillatory relationships 

clear, clc, clf
X = randn(1,3000);
% Y = randn(1,3000);

% Shift positions of elements circularly
Y = circshift(X,124);

h1 = subplot(211);
plot(X)
xlabel('Time (ms)')

h2 = subplot(212);
plot(Y)
xlabel('Time (ms)')

linkaxes([h1 h2],'x')

% remember that the reference signal is the second
% input of the xcorr function  

[CCG, lags] = xcorr(Y,X);

plot(lags,CCG)

xlabel('Lag (ms)')
ylabel('CCG')

% ylim([-3000 3000])


%% another exemple

clear, clc, clf

Stimulus = zeros(1,1000);

Istim = randi(975,[1,20]);

Stimulus(Istim) = 1;

Y = randi(40,[1,1000]);

% Y(Istim+20) = 0*Y(Istim+20);
% Y(Istim+21) = 0*Y(Istim+21);
% Y(Istim+22) = 0*Y(Istim+22);
% Y(Istim+23) = 0*Y(Istim+23);
% Y(Istim+24) = 0*Y(Istim+24);

Y(Istim+20) = 2*Y(Istim+20);
Y(Istim+21) = 2*Y(Istim+21);
Y(Istim+22) = 2*Y(Istim+22);
Y(Istim+23) = 2*Y(Istim+23);
Y(Istim+24) = 2*Y(Istim+24);

h1 = subplot(211);
plot(Y)
xlabel('Time (ms)')

h2 = subplot(212);
plot(Stimulus)
xlabel('Time (ms)')

linkaxes([h1 h2],'x')

%%
[CCG, lags] = xcorr(Y,Stimulus,100);

clf
plot(lags,CCG)
xlabel('Lag (ms)')
ylabel('CCG')

%% Oscilatory Case

clear

RandomSpikes = randi(1000,[1,100]);
SpikeTime(RandomSpikes) = 1;
SpikeTime(20:20:1000) = 1;
subplot(2,1,1)
plot((1:1000)*1/1000,SpikeTime)

[ACG, lags] = xcorr(SpikeTime,200);
ACG(lags==0) = 0;
subplot(2,1,2)
plot(lags,ACG)

xlabel('Lag (ms)')
ylabel('ACG')

%% Oscilatori Cases from two Signals
% i.e., cross-correlation

clear

srate = 1000;
dt = 1/srate;
t = dt:dt:3;

LFP1 = 1*sin(2*pi*4*t)+0*sin(2*pi*30*t)+0*randn(size(t));
LFP2 = 0*sin(2*pi*4*t+pi)+1*sin(2*pi*77*t)+0*randn(size(t));

subplot(211)
plot(t,LFP1)
hold on
plot(t,LFP2)
hold off
% xlim([1 1.5])

subplot(212)
[CCG, lags] = xcorr(LFP2,LFP1);
plot(lags*dt,CCG)

xlabel('Time (s)')
ylabel('CCG')

% ylim([-1500 1500])

%% applications for real cases

clear

load('LFP_HG_HFO.mat')

srate = 1000;
dt = 1/srate;
t = dt*(1:length(lfpHG));


% [ACG lags] = xcorr(lfpHG,1000);
[CCG, lags] = xcorr(lfpHG,lfpHFO,1000);

subplot(211)
plot(t,lfpHG)
xlim([3 5])

subplot(212)
% plot(lags*dt,ACG)
plot(lags*dt,CCG)


%% normalizing the CCG

% if you want the CCG to express the Pearson
% correlation coefficient, just compute the CCG
% between the z-scores of the time series

clear

srate = 1000;
dt = 1/srate;
t = dt:dt:3;

LFP1 = 1*sin(2*pi*4*t)+1*sin(2*pi*30*t)+0*randn(size(t));
LFP2 = 1*sin(2*pi*4*t+0*pi)+1*sin(2*pi*77*t)+0*randn(size(t));

normLFP1 = zscore(LFP1);
normLFP2 = zscore(LFP2);

clf

% [CCG lags] = xcorr(LFP2,LFP1)

[CCG, lags] = xcorr(normLFP2,normLFP1);
CCG_pearson = CCG/(length(LFP1)-1);
plot(lags*dt,CCG_pearson)

xlabel('Time (s)')
ylabel('Pearson correlation')

% [CCG lags] = xcorr(normLFP2/sqrt(length(normLFP2)),...
%     normLFP1/sqrt(length(normLFP1)))
% 
% plot(lags*dt,CCG)
% xlabel('Time (s)')
% ylabel('Pearson correlation')

ylim([-1 1])

%%

[r, p] = corr(LFP1(66:end)',LFP2(1:end-65)');

%% Matlab built in normalization for CCG

[CCGB, lags] = xcorr(LFP2,LFP1,'coeff');

hold on

plot(lags*dt,CCGB)

hold off


%% Relationship between Power and ACG

LFP = sin(2*pi*8*t);

[ACG, lags] = xcorr(LFP);

subplot(211)
plot(lags*dt,ACG)

% the Fourier transform at frequency f is
% equivalent to power 

Power = zeros(1, 20);
for f=1:20
   Power(f) = sum(ACG.*exp(-1i*2*pi*f*(lags*dt))); 
end
    
subplot(212)
plot(1:20,real(Power),'bo-')
hold on
plot(1:20,imag(Power),'ro-')
hold off

























