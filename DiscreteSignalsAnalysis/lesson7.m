%% Lesson 7 = 2021/08/24; STA - Spike-Triggered Average

clear, clc, clf
load('LFPBuz.mat')

%% playing with sines

% setup
clear, clf
f=8;
srate=1000;
dt=1/srate;
t=dt:dt:10;

% create fake LFP - sine qith frequency f
LFP = sin(2*pi*f*t);

% set index of periodic spike signal with
% consistent relation with the sine wave in the 8
% Hz LFP
Idxmodulado = 125:125:10000;
% Plus some random spike indices
Idxruido = randi(10000,[1,200]);
Idx = [Idxmodulado, Idxruido];

% assign 1 to every point were there is a spike
spkind(Idx)=1;

subplot(211)
    plot(t,LFP);         hold on
    plot(t(Idx),0,'ko'); hold off
    xlim([0 2])
    xlabel('Time (s)')
    ylabel('mv')
    title('Spike activity over LFP')

% constructing the STA
% I need to be recomputed because some of the
% random spikes can have coincided with the
% modulated ones
I = find(spkind==1);
Nspike = length(I);

winl = 0.6*srate; % num of points in win
clear STA
STA2 = zeros(1,winl+1);
count = 0;
for nspike=1:Nspike
    % if the given spike is not so much on the
    % edge of the signal that there isn't an
    % entire window around it
if (I(nspike)>winl/2 & I(nspike)<length(LFP) - winl/2)
    count=count+1;
    % take all points around spike that fit to win
    winidx = I(nspike) - round(winl/2):I(nspike)+round(winl/2);
    
    % for visual inspection: plot
    % hold on
    % plot(t(winidx),LFP(winidx),'y')
    % hold off

    %%% two ways of doing it
    % store all win values (more memory intensive)
    STA(count,:) = LFP(winidx);
    % only sum all points cause we are only interested in the mean
    STA2 = STA2 + LFP(winidx); % sum, not concatanation
end
end

% complete/take the mean
STA2 = STA2/count;
STA = mean(STA);

subplot(212)
    plot(((1:length(STA))-winl/2 )*dt,STA)
    hold on
    plot([0 0],[-max(STA) max(STA)],'k--')
    hold off
    xlim([-winl/2 winl/2]*dt)
    ylim([-1 1])
    xlabel('Time (s)')
    ylabel('mv')
    title('STA')

%% Showing how the number of spikes, if independent of LFP phase, can influence STA

% see de changes by making n bigger os smaller
STA = zeros(1,srate);
subplot(3,1,[1 2])
    for n=1:30
        Sinal=sin(2*pi*8*t(1:srate)+2*pi*rand);

        plot(t(1:srate)-0.5,Sinal+n,'k')
        hold on

       STA = STA+Sinal;
    end
    STA = STA/n;
    hold off
    xlabel('Time (s)')
    ylabel('LFP #')
    title('STA')
    xlim([0 1]-0.5)

subplot(313)
    plot(t(1:srate)-0.5,STA)
    xlabel('Time (s)')
    ylabel('mv')
    title('STA')
% ylim([-1 1])


%% Programing Surrogates

% set up
clear, clc, clf
srate=1000;
dt = 1/srate;
t = dt:dt:10;
f = 8; % in Hz
LFP = sin(2*pi*f*t);

Idxmodulado = 125:125:10000;
Idxruido = randi(10000,[1,20]);
Idx = [Idxmodulado,Idxruido];
spkind(Idx)=1;

subplot(311)
    plot(t,LFP); hold on
    plot(t(Idx),LFP(Idx),'ko'); hold off
    xlim([0 10])
    xlabel('Time (s)')
    ylabel('mv')
    title('Spike activity over LFP')

% constructing the STA (same as before)
I = find(spkind==1);
Nspike = length(I);
winl = 0.6*srate;
STA = zeros(1,winl+1);
count = 0;
for nspike=1:Nspike
if I(nspike) > winl/2 & I(nspike)< length(LFP) - winl/2
    count = count+1;
    winidx = I(nspike)- round(winl/2):I(nspike)+round(winl/2);
    STA = STA + LFP(winidx);
end
end
STA = STA/count;

subplot(312)
    plot(((1:length(STA))-winl/2 )*dt,STA); hold on
    plot([0 0],[-max(STA) max(STA)],'k--');      hold off
    xlim([-winl/2 winl/2]*dt)
    ylim([-1 1])
    xlabel('Time (s)')
    ylabel('mv')
    title('STA')

STAreal = STA;

% computing the STA chance
% get 100 surrogates STA each avg of Nspikes random spikes
clear STAsurr
for nsurr = 1:100
    STA = zeros(1,winl+1);
    count = 0;   
    for nspike=1:Nspike
        Isurr = randi(length(LFP));
        if Isurr > winl/2 & Isurr< length(LFP) - winl/2
            count = count+1;
            winidx = Isurr- round(winl/2):Isurr+round(winl/2);
            STA = STA + LFP(winidx);  
        end
    end   
    STAsurr(nsurr,:) =   STA/count;
end

subplot(313)
    plot(((1:length(STAreal))-winl/2 )*dt,STAsurr,'color',[1 1 1]/1.3); hold on
    plot(((1:length(STAreal))-winl/2 )*dt,STAreal,'b-','linew',3)
    plot(((1:length(STAreal))-winl/2 )*dt,mean(STAsurr),'color',[0 0 0])
    plot(((1:length(STAreal))-winl/2 )*dt,mean(STAsurr)+3*std(STAsurr),'k--')
    plot(((1:length(STAreal))-winl/2 )*dt,mean(STAsurr)-3*std(STAsurr),'k--')
    plot([0 0],[-max(STAreal) max(STAreal)],'k--'); hold off
    xlim([-winl/2 winl/2]*dt)
    xlabel('Time (s)')
    ylabel('mv')
    title('STA')
% ylim([-1 1])


%% Applying to real data 

% set up
clear, clf, clc
load('SpkBuz.mat')
load('LFPBuz.mat')

%create time vector
dt = 1/srate; % sampling period, in seconds (s)
t = dt*(1:length(LFP)); % in s

neuron = 7;
spktimes = Raster{neuron}; % in s

I = round(spktimes*srate); % round ms to indexes
Nspike = length(I);

subplot(211)
    plot(t,LFP)
    hold on
    % plot(spktimes,0,'ko')
    plot(t(I),LFP(I),'ko','markerf','k')
    hold off
    xlabel('Time (s)')
    ylabel('mv')
    title('Spike activity over LFP')

% constructing the STA (same as before)
winl = 2*srate;
STA = zeros(1,winl+1);
count = 0;
for nspike=1:Nspike
if I(nspike) > winl/2 & I(nspike)< length(LFP) - winl/2
    count = count+1;
    winidx = I(nspike)- round(winl/2):I(nspike)+round(winl/2);
    STA = STA + LFP(winidx)';
end
end
STA = STA/count;
STAreal = STA;

% %
subplot(212)
    plot(((1:length(STA))-winl/2 )*dt,STA); hold on
    plot([0 0],[min(STA) max(STA)],'k--'); hold off
    xlim([-winl/2 winl/2]*dt)
    xlabel('Time (s)')
    ylabel('mv')
    title('STA')

% computing the STA chance
% get 100 surrogates STA each avg of Nspikes random spikes
clear STAsurr
for nsurr = 1:100
    STA = zeros(1,winl+1);
    count = 0;   
    for nspike=1:Nspike
        Isurr = randi(length(LFP));
        if Isurr > winl/2 & Isurr< length(LFP) - winl/2
            count = count+1;
            winidx = Isurr- round(winl/2):Isurr+round(winl/2);
            STA = STA + LFP(winidx)';  
        end
    end
    STAsurr(nsurr,:) = STA/count;
end

% %
subplot(212)
    plot(((1:length(STAreal))-winl/2 )*dt,STAsurr,'color',[1 1 1]/1.3); hold on
    plot(((1:length(STAreal))-winl/2 )*dt,STAreal,'b-','linew',3)
    plot(((1:length(STAreal))-winl/2 )*dt,mean(STAsurr),'color',[0 0 0])
    plot(((1:length(STAreal))-winl/2 )*dt,mean(STAsurr)+3*std(STAsurr),'k--')
    plot(((1:length(STAreal))-winl/2 )*dt,mean(STAsurr)-3*std(STAsurr),'k--')
    plot([0 0],[min(STAreal) max(STAreal)],'k--'); hold off
    xlim([-winl/2 winl/2]*dt)
    xlabel('Time (s)')
    ylabel('mv')
    title('STA')
    % ylim([-1 1])




















