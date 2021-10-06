%% Lesson 9 = 06/04/2021; EEGFILT - Filters

% filtring using eegfilt function

clear
clc
clf

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

LFP = sin(2*pi*10*t)+sin(2*pi*30*t)+sin(2*pi*50*t);

% eegfilt function calls filtfilt which applies
% the filter twice (a reverse to correct phase
% distortion) 

% ordem default of the filter is given by:

% default 3*fix(srate/locutoff)

% Low-pass filter (filtro passa baixa)
% filtrado = eegfilt(LFP,srate,0,16);

% High-pass filter (filtro passa alta)
% filtrado = eegfilt(LFP,srate,35,0);

% Band-pass filter (filtro passa banda)
filtrado = eegfilt(LFP,srate,45,55);

% changing filter order
% ordem = 66
% filtrado = eegfilt(LFP,srate,45,55,0,ordem);

% Notch filter (filtro rejeita banda)
% filtrado = eegfilt(LFP,srate,20,40,0,[],1);


[Pxx , ~] = pwelch(LFP,length(LFP),[],2^16,srate);
[PxxFilt, F] = pwelch(filtrado,length(LFP),[],2^16,srate);

figure(2)

subplot(311)
plot(t,LFP)
hold on
plot(t,filtrado,'r-')
hold off

subplot(312)
plot(F,Pxx)
hold on
plot(F,PxxFilt,'r-')
hold off
xlim([0 70])

subplot(313)
plot(F,PxxFilt,'r-')
xlim([0 70])


%% No filter is Perfect

% clf

figure(1)
% clf

LowFreqCutoff=55;
HighFreqCutoff=65;

whitenoise = randn(1,1000000);
 
% filtered = eegfilt(whitenoise,srate,...
%     LowFreqCutoff,HighFreqCutoff);
% 
% ordem = 100
% filtered = eegfilt(whitenoise,srate,...
%     LowFreqCutoff,HighFreqCutoff,0,ordem);

filtered = eegfilt(whitenoise,srate,...
    LowFreqCutoff,HighFreqCutoff,0,[],1); % change the last by 0

[PxxW, ~] = pwelch(whitenoise,srate,[],2^16,srate);
[PxxF, F] = pwelch(filtered,srate,[],2^16,srate);

subplot(211)
plot(F,PxxW,'k-')
hold on
plot(F,PxxF)
plot([HighFreqCutoff HighFreqCutoff],[0 max(PxxW)*1.2],'k-')
plot([LowFreqCutoff LowFreqCutoff],[0 max(PxxW)*1.2],'k-')

hold off

% xlim([0 70])

xlim([0 100])

%% understanding the fifth imput of eegfilt

% epochframes = frames per epoch 
%(filter each epoch separately {def/0: data is 1 epoch}

t = dt:dt:1;

LFP1 = sin(2*pi*10*t);
LFP2 = sin(2*pi*10*t+pi);
LFP3 = sin(2*pi*10*t+pi/2);
LFP4 = sin(2*pi*10*t+pi/3);

LFP = [LFP1 LFP2 LFP3 LFP4];

filtradoTrials = eegfilt(LFP,srate,0,15,1000);
filtradoSingle = eegfilt(LFP,srate,0,15);

clf

plot(LFP,'k-')
hold on
plot(filtradoSingle)
plot(filtradoTrials,'r-')
hold off

%% Effect of order size again (no leakage temporal)

srate = 1000;
dt = 1/srate;
t = dt:dt:4;

lfp = sin(2*pi*10*t);
lfp(1:2000)=0;
lfp(3000:end)=0;

ordem = 500;
filtrado = eegfilt(lfp,srate,0,15,0,ordem);

plot(t,lfp)
hold on
plot(t,filtrado,'r-','linew',2)
hold off


xlim([1.5 3.5])

%% No filter is perfect 2


srate = 1000;
dt = 1/srate;
t = dt:dt:4;

lfp = sin(2*pi*10*t);
lfp = lfp + 0.0*randn(size(t));

filtrado = eegfilt(lfp,srate,15,0);

[Pxx , ~] = pwelch(lfp,4*srate,[],2^16,srate);
[PxxF, F] = pwelch(filtrado,4*srate,[],2^16,srate);


subplot(221)
plot(t,lfp)
hold on
plot(t,filtrado)
hold off
xlim([1 1.5])

subplot(222)
plot(F,Pxx)
hold on
plot(F,PxxF)
hold off
xlim([0 20])

subplot(223)
plot(t,filtrado,'r-')
xlim([1 1.5])

subplot(224)
plot(F,PxxF,'r-')
xlim([0 20])

%% The eegfilt accept matrices as inputs

% in this case it filter per line

% Channel vs timepoints (that is, every line is a
% channel)

LFPAll(1,:) = sin(2*pi*5*t);
LFPAll(2,:) = sin(2*pi*15*t);
LFPAll(3,:) = sin(2*pi*30*t);

filtered = eegfilt(LFPAll,srate,10,25,0,[],1);


h1 = subplot(311);
plot(t,LFPAll(1,:))
hold on
plot(t,filtered(1,:))
hold off
ylim([-3 3])

h2 = subplot(312);
plot(t,LFPAll(2,:))
hold on
plot(t,filtered(2,:))
hold off
ylim([-3 3])

h3 = subplot(313);
plot(t,LFPAll(3,:))
hold on
plot(t,filtered(3,:))
hold off
ylim([-3 3])

linkaxes([h1 h2 h3])

% linkaxes([h1 h2 h3],'x')

%%


















