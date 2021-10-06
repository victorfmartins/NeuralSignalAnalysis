%% Lesson 4 = 2021/08/03; Rastergrams and PSTH (peri-stimulus time histogram)

clear, clc
load('SpikeGustatoryCortex.mat')

% dimensions (of data spiketimes just loaded)
% LGC = left gustatory cortex; RGC = right GC
% 1 - Hemisphere (LGC: 1, RGC: 2)
% 2 - Mouse number (mouse: 1, 2 e 3)
% 3 - Session (1-4, experiment days)
% 4 - Taste (1: NaCl, 2: Citric Acid, 3: Quinine (bitter), 4: Sucrose)
% 5 - Unit Number (max of 8 neurons)
% 6 - Trial Number (max of 10 trials)

%% Single spike train rastergram

% setup - data we want
clear, clc, clf
load('SpikeGustatoryCortex.mat')
hemi = 1;
rat = 1;
session = 1;
taste = 1;
unit = 1;
trial = 1;

% taking one second out from basal activity (mark input time)
spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 

subplot(111)
    plot(spks,ones(size(spks)),'ko','markerfacecolor','c')
    hold on
    plot([0 0],[0 2],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    set(gca,'ytick',1,'fontsize',14)
    xlabel('Time (s)')
    % ylabel('LGC, rat 1, day 1, Taste NaCl, Neuron 1, Trial #')
    ylabel('Trial #')
    title('Rastergram')


%% Loop over trials - Rastergram

% setup - data we want
clf
hemi = 1;
rat = 1;
session = 1;
taste = 1;
unit = 1;

subplot(111)
    for trial = 1:10
        spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 
        plot(spks,trial*ones(size(spks)),'ko','markerfacecolor','c')
        hold on
    end
    plot([0 0],[0 trial+1],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    ylim([0 10])
    set(gca,'ytick',1:trial,'fontsize',14)
    xlabel('Time (s)')
    ylabel('Trial #')
    title('Rastergram')


%% Loop over tastes and trials- Rastergram

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

%% Discrete PSTH

% setup - data we want
clf
hemi = 1;
rat = 1;
session = 1;
unit = 1;

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
    set(gca,'fontsize',14)
    title('Rastergram')

    
%%% Discrete PSTH
winlength = 0.1;
clear PSTH T
nwintotal = 3.5/winlength;
for nwin = 1:nwintotal
    LimE = -1 + (nwin-1)*winlength;
    LimD = -1+winlength + (nwin-1)*winlength; 
    PSTH(nwin) = sum(SPKSALL>LimE & SPKSALL<=LimD);
    T(nwin) = (LimE+LimD)/2;
end

subplot(3,1,3)
    % transforming count in Firing Rate (in Hz)
    bar(T,PSTH/(winlength*count),'k')
    set(gca,'fontsize',14)
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')


%% Continuous PSTH

clear PSTH T
%%% This makes the continuous aspect. it has to be step<=winlength. 
%%% If step=winlength then this plot is discrete rather then continuos
%%% If step<windlength then the continuous aspect starts to be seem
%%% try step = 0.01 to see continuous aspect. 
%%% try step<0.01 to see that nothing happens
step = 0.1;
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
    ylabel('Firing Rate (Hz)')
    set(gca,'fontsize',14)
    
    
%% Discrete PSTH without loop rotine

% get all spikes without loop
SPKSALL2 = [spiketimes{hemi,rat,session,:,unit,:}]-1;
bincenters = -0.95:winlength:2.45;
PSTH2 = hist(SPKSALL2,bincenters);

subplot(3,1,3)
    bar(T,PSTH,'k') %%% notice spike count here, not firing rate
    hold on
    % bar(T,PSTH2,'g')
    hold off
    set(gca,'fontsize',14)
    xlabel('Time (s)')
    ylabel('spike counts')
    
    
%% (Discrete) PSTH with standard deviation from trials bar

clf
clear FR FR2
clear PSTH T T2
winlength = 0.1;
step = 0.1;
Nwin = (3.5-winlength)/step+1;

subplot(3,1,[1 2])
    count = 0;
    for taste = 1:4
    for trial = 1:10
        count = count+1;
        spks = spiketimes{hemi,rat,session,taste,unit,trial}-1; 
        plot(spks,(trial+(taste-1)*10)*ones(size(spks)),...
        'ko','markerfacecolor',cor{taste})
        hold on
        
        % discrete case
        [spikecount T] = hist(spks,-0.95:winlength:2.45); 
        FR(count,:) = spikecount/winlength;
        
        % continuous case
        for nwin = 1:Nwin
           LimE = -1 + (nwin-1)*step;
           LimD = LimE + winlength;
           PSTH(nwin) = sum(spks>LimE & spks<=LimD);
           T2(nwin) = (LimE+LimD)/2;   
        end
        FR2(count,:) = PSTH/winlength;
    end
    end
    plot([0 0],[0 count+1],'r--','linew',2)
    hold off
    xlim([-1 2.5])
    ylim([0 count])
    ylabel('Trial #')
    set(gca,'fontsize',14)
    title('Rastergram')

subplot(3,1,3)
    M = mean(FR);
    SD = std(FR);
    SEM = SD/sqrt(count);
    errorbar(T,M,SEM,'k-sq','markerfacecolor','k')
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    set(gca,'fontsize',14)

%% (Continuous) PSTH with standard deviation from trials bar

M2 = mean(FR2);
SD2 = std(FR2);
SEM2 = SD2/sqrt(count);

subplot(313)
    plot(T2,mean(FR2),'k','linew',2)
    hold on
    plot(T2,mean(FR2)+SEM2,'k','linew',1)
    plot(T2,mean(FR2)-SEM2,'k','linew',1)
    hold off

    A = mean(FR2)+SEM2;
    B = mean(FR2)-SEM2;

    fill([T2 T2(end:-1:1)],[A B(end:-1:1)],'c')
    hold on
    plot(T2,mean(FR2),'k','linew',2)
    hold off
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    set(gca,'fontsize',14)










