%% Lesson 1 = 2021/07/22; Spike detection
% Spike detenction is a holl field on its on. If
% one venture in to details one may get a PhD on it

% Paper recomendation: Spike sorting with
% Gaussian mixture model - Bryan de Souza

% For this class we used tetrodo.mat data. This is
% a 4 very close to each other eletrode bondle.

clear
clc
clf

I = imread('Tetrodo.jpg');
imshow(I)
clc
% Figure legend:
% TetrodLegend = ['Multisite electrodes (a wire tetrode, for example)' ...
%     newline 'can estimate the position of the recorded neurons by' ...
%     newline 'triangulation. Distance of the visible electrode tips' ...
%     newline 'from a single pyramidal cell (triangles) is indicated' ...
%     newline 'by arrows. The spike amplitude of neurons (>60 ?V)' ...
%     newline 'within the gray cylinder (50 ?m radius), containing' ...
%     newline '?100 neurons, is large enough for separation by' ...
%     newline 'currently available clustering methods. Although the' ...
%     newline 'extracellularly recorded spike amplitude decreases'
%     newline 'rapidly with distance, neurons within a radius of' ...
%     newline '140 ?m, containing ?1,000 neurons in the rat cortex,' ...
%     newline 'can be detected. Improved recording and clustering' ...
%     newline 'methods are therefore expected to record from larger' ...
%     newline 'number of neurons in the future. (Data are derived' ...
%     newline 'from simultaneous extracellular and intracellular'...
%     newline 'recordings from the same pyramidal cells.)']
load('tetrodo.mat','data')

srate = 30000; % sampling rate, in Hz (1/s)
dt = 1/srate; % sampling period, in seconds (s)
t = (1:length(data))*dt;

%% Show data from one channel

ch = 1;
rawsignal = data(:,ch);
plot(t,rawsignal)
ylabel('mV')
    xlabel('Time (s)')
xlim([0 60])


%% Signal TV

stepsize = 0.01;

for nstep = 0:100
xlim([0 0.2]+nstep*stepsize)    
    ylim([-10 10])
    pause(0.001)
end

%% Geting MUA - Multi-Unit Activity
% Frequencies cann't be too low here

lowf = 500;
highf = 4000;
    
MUA = eegfilt(rawsignal',srate,lowf,highf);

%% Plotando MUA e LFP
% The spikes have les than 1mV and can be seen in
% the LFP with the proper zoom in.

h1 = subplot(2,1,1);
plot(t,rawsignal,'k')
title('Raw Signal')
ylabel('mV')
xlabel('Time (s)')
xlim([0 60])
h2 = subplot(2,1,2);

plot(t,MUA,'k')
title('MUA')
ylabel('mV')
xlabel('Time (s)')
xlim([0 60])

linkaxes([h1 h2],'x')
    
%% Computing the MUA of all eletrodes

MUA1 = eegfilt(data(:,1)',srate,lowf,highf);
MUA2 = eegfilt(data(:,2)',srate,lowf,highf);
MUA3 = eegfilt(data(:,3)',srate,lowf,highf);
MUA4 = eegfilt(data(:,4)',srate,lowf,highf);

%% Ploting the MUAs with a vertical displacement
% one can see spike wave shapes by zooming on one
% spike to 5milisecond frame

% each eletrode can yield a differente wave form
% for a given spike

clf

plot(t,MUA1,'k')
hold on
plot(t,MUA2-1,'b')
plot(t,MUA3-2,'r')
plot(t,MUA4-3,'g')
hold off
axis tight

%% Starting over with different data
% one can seen waveshape (signature) of diferent
% neurons with a proper zooming to 5miliseconds

clear 

% j? ? MUA:
load('rawdata3.mat')

srate = 20000; % in Hz
dt = 1/srate;
t = dt*(1:length(data)); % in seconds

subplot(111)
    plot(t,data)
    xlabel('Time (s)')
    ylabel('mV')
    title('MUA')

%% Gathering up pointing peaks
% Define voltage treshold to classify a point as
% belonging to a spike
% Stain all points above treshold with red
% Get the upper point in each spike
% WARNING: this rotine may not catch downward
% spikes

% arbitrary treshold definition (anything well
% above base activity and below usual max spike
% voltage is well suited. Too low treshold leads
% to spike classification of non spike flutuation.
% Too high treshold leads to missis in spike
% detection. Finding all adequate boundaries is a
% matter of a Phd thesis.
% Here detrending is necessary if one have trended data.
threshold = mean(data)+2*std(data);

plot(t,data,'b.-')
hold on
% plot treshold line
plot([0 t(end)],[threshold threshold],'y-','linew',2)
hold off

% identify pointss to go red
I = find(data>threshold);

hold on
% red identified points
plot(t(I),data(I),'r.')
hold off
% a peak point has a higger value than both its
% neiboors
nspikes = 0;
clear Ispikes
for j=2:length(I)-1
    if data(I(j))>data(I(j-1))
        if data(I(j))>data(I(j+1))
            nspikes = nspikes+1;
            Ispikes(nspikes)=I(j);
        end
    end
end
hold on
% black triangle peak staining
plot(t(Ispikes),data(Ispikes),'kv','markerf','k','markersize',10)
hold off
    xlabel('Time (s)')
    ylabel('mV')
    title('MUA')

%% Findpeaks with built in matlab function

[temp,Ispikes2]= findpeaks(data,'MinPeakHeight',threshold);
hold on
plot(t(Ispikes2),data(Ispikes2),'ko','markerf','r','markersize',10)
hold off

%% Staining waveshape
% staining neiboors of a peak point

hold on
for j=1:length(Ispikes2)
    % 40 neiboors to each side
    plot(t(Ispikes2(j)-40:Ispikes2(j)+40),...
        data(Ispikes2(j)-40:Ispikes2(j)+40),'g-')
end
hold off

%% Center the wavepeak and plot all waveshapes
% gather neiboors of a peak point
clear waveshape

% pré alocation of memori for claimed speed
% % waveshape = zeros(length(Ispikes2),81);
tic
for j=1:length(Ispikes2)
    % 20 neiboors to each side
    waveshape(j,:) = data(Ispikes2(j)-20:Ispikes2(j)+20);
end
toc
plot(-20:20,waveshape','k')

%%% Next move is to use all the neiboors points of
%%% each spike to feed a spike sorting algoritm.
%%% See next class.

























