%% 10/04/2021 - Exam Resolution
%/ ************************************************************************** #
%                                                                             #
%                                              ::::::  ::::::    :::    :::   #
%  Exam 1 Resolution                        :+:     :+:         :+:+:  :+:+:  #
%                                          +:+     +:+         +:+ :+ :+ +:+  #
%  By: vfranco- <victorf.martins@usp.br>  +#+     +#+         +#+   +:+ +#+   #
%                                        +#+     +#+         #+#       #+#    #
%  Created: 2021/04/10 13:11              #+#     #+#       #+#       #+#     #
%  Updated: 2021/04/16 00:26                ######  ###### ###       ###.usp  #
%                                                                             #
% *************************************************************************** #

%% Question 1
% 1) Plot all channels simultaneously (no overlap).

clc % clean comm win
clf % clean fig
clear % clean mem
load('LFPprobe.mat')

dt = 1/srate;  % in s
Tmax = length(LFPprobe)/srate;  % in s
t = dt:dt:Tmax;  % in s

[macs,~] = max(abs(LFPprobe),[],2);
m = max(macs);
clear macs

altura = 0:100:1600;

% mode C for collor, mode Z for zoom and black
mode = 'Z';

figure(1)
if mode == 'Z'
    for i = 1:16
    plot(t,LFPprobe(i,:)+i*1000,'k-','LineWidth',1.2)
    hold on
    end
    hold off
    set(gca,'ytick',(0:16)*1000,'yticklabel',altura)
    xlim([96, 98])
    axis ij
    xlabel('Time (s)')
    ylabel('Profundidade (\mum)')
else
    for i = 1:16
        plot(t,LFPprobe(i,:)+i*m)
        hold on
    end
    hold off
    set(gca,'ytick',(0:16)*m,'yticklabel',altura)
    axis ij
    xlabel('Time (s)')
    ylabel('Profundidade (\mum)')
    xlim([0 600])
end

%% Question 2
% 2) Vary the period shown using 2 s sliding
% windows (“LFPTv”). 

clc, clf, clear
load('LFPprobe.mat')

dt = 1/srate; % in s
Tmax = length(LFPprobe)/srate; % 600 s (10 min)
t = dt:dt:Tmax; % in s

% [macs,~] = max(abs(LFPprobe),[],2);
% m = max(macs);
% clear macs

altura = 0:100:1600;

for i = 1:16
    plot(t,LFPprobe(i,:)+i*3000,'k-','LineWidth',1)
    hold on
end

hold off
set(gca,'ytick',(0:16)*3000,'yticklabel',altura)
xlim([0 1])
axis ij
xlabel('Time (s)')
ylabel('Profundidade (\mum)')


stepsize = 0.2;
for nstep = 0:100
    xlim([0 1]+nstep*stepsize)
    pause(0.001)
end


%% Question 3
% 3) Compute the PSD of each channel, both in
% absolute values and normalized by the total
% power (% power).  

clc, clf, clear, close
load('LFPprobe.mat')

dt = 1/srate; % in s
Tmax = length(LFPprobe)/srate; % 600 s (10 min)
t = dt:dt:Tmax; % in s

win = 1*srate;
overlap = [];
nfft = 10000;

for i  = 1:16
[Pxx(:,i), F(:,i)] = pwelch(LFPprobe(i,:),win,overlap,nfft,srate);

% Percentage power
NormPxx(:,i) = Pxx(:,i)/sum(Pxx(:,i));
end

% EXTRA
% figure(3)
% window_length = 2*srate;
% step_size = 0.2*window_length; % (ie., 80% overlap)
% 
% Nwindows = (length(LFPprobe(1,:))-window_length)/step_size+1;
% 
% clear T TFD
% 
% for nwin = 1:Nwindows
%  winidx = (1:window_length) + (nwin-1)*step_size;
% [Pxx F] = pwelch(LFPprobe(1,(winidx)),window_length,[],nfft,srate);
% T(nwin) = t(winidx(window_length/2));
% TFD(nwin,:) = Pxx; 
%     
% end
% % %
% 
% imagesc(T,F,TFD')
% axis xy
% ylim([0 30])
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% 
% % colorbar
%% Question 4
% 4) Plot four graphs in distinct subplots, where
% each graph shows the PSDs of all channels,
% according to the following variations: [absolute
% or normalized] vs [linear or logarithmic scale].   

figure(4)
subplot(2,2,1)
plot(F,Pxx)
xlabel('Frequency (Hz)')
ylabel('Power mV^2/Hz')
xlim([0 20])

subplot(2,2,2)
plot(F,NormPxx)
xlabel('Frequency (Hz)')
ylabel('Normalized Power')
xlim([0 20])

subplot(2,2,3)
plot(F,10*log10(Pxx))
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0 100])
ylim([-10 50])

subplot(2,2,4)
plot(F,10*log10(NormPxx))
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0 100])
ylim([-70 -10])


%% Question 5
% 5) Compute the average power (absolute and
% normalized) in the theta band (5-10 Hz) for each
% channel.  

ITheta = find(F(:,1) > 5 & F(:,1)<10);

PxxThetaPower = mean(Pxx(ITheta,:));
% Percentage power
NormThetaPower = mean(NormPxx(ITheta,:));


%% Question 6
% 6) Plot the values obtained in step 5 as a line
% graph, indicating on the Y axis the location
% (“height”) of the electrode (consider 100 to
% 1600 mum) and on the X axis the average power
% value. Plot one for absolute values and one for
% normalized values.      

altura = 100:100:1600;

figure (6)
subplot(1,2,1)
plot(PxxThetaPower,altura,'k^-','MarkerFaceColor','k')
axis ij
xlabel('\Theta Power')
ylabel('Profundidade (\mum)')

subplot(1,2,2)
plot(NormThetaPower,altura,'k^-','MarkerFaceColor','k')
axis ij
xlabel('Norm \Theta Power')
ylabel('Profundidade (\mum)')


%% Question 7
% 7) Repeat steps 5 and 6 for the gamma (30 to 100
% Hz) and HFO (120-160 Hz) bands. 


%%%% STRETCH THE GRAPH HORIZONTALLY TO SEE BETTER %%%%


%%%%%% gama %%%%%%
IGama = find(F(:,1) > 30 & F(:,1)<100);

% Percentage power
PxxGamaPower = mean(Pxx(IGama,:));
NormGamaPower = mean(NormPxx(IGama,:));

figure(7)
subplot(1,4,1)
plot(PxxGamaPower,altura,'k^-','MarkerFaceColor','k')
axis ij
xlabel('Gama Power')
ylabel('Profundidade (\mum)')

subplot(1,4,2)
plot(NormGamaPower,altura,'k^-','MarkerFaceColor','k')
axis ij
xlabel('Norm Gama Power')
ylabel('Profundidade (\mum)')

%%%%%%% HFO %%%%%%
IHFO = find(F(:,1) > 120 & F(:,1)<160);

% Percentage power
PxxHFOPower = mean(Pxx(IHFO,:));
NormHFOPower = mean(NormPxx(IHFO,:));

subplot(1,4,3)
plot(PxxHFOPower,altura,'k^-','MarkerFaceColor','k')
axis ij
xlabel('HFO Power')
ylabel('Profundidade (\mum)')

subplot(1,4,4)
plot(NormHFOPower,altura,'k^-','MarkerFaceColor','k')
axis ij
xlabel('Norm HFO Power')
ylabel('Profundidade (\mum)')

%% Question 8
% 8) Compute and plot the TFD (spectrogram) for
% channel 1 and channel 16 in distinct subplots. 

% factor that determines the trade off frqXtime
window_length = 2*srate;

step_size = 0.2*window_length; % (ie., 80% overlap)
Nwindows = (length(LFPprobe(1,:))-window_length)/step_size+1;


%%%% Normalized
figure(8)
subplot(211)
    clear T TFD

    for nwin = 1:Nwindows
     winidx = (1:window_length) + (nwin-1)*step_size;
    [Pxx, F] = pwelch(LFPprobe(1,(winidx)),...
                window_length,[],nfft,srate);
    T(nwin) = t(winidx(window_length/2));
    TFD1(nwin,:) = Pxx/sum(Pxx);    
    end

    imagesc(T,F,TFD1')
    axis xy
    ylim([0 30])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Channel 1')
    colorbar % can uncomment to see color scales

subplot(212)
    clear T TFD

    for nwin = 1:Nwindows
     winidx = (1:window_length) + (nwin-1)*step_size;
    [Pxx, F] = pwelch(LFPprobe(16,(winidx)),...
                window_length,[],nfft,srate);
    T(nwin) = t(winidx(window_length/2));
    TFD16(nwin,:) = Pxx/sum(Pxx);    
    end

    imagesc(T,F,TFD16')
    axis xy
    ylim([0 30])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Channel 16')
    colorbar % can uncomment to see color scales

%% Question 9
% 9) Use the result in 8 to compute and plot the
% time series of the average power in theta, gamma
% and HFO in the two channels (show each band in
% distinct subplots).    

figure(9)
subplot(311)
    ThetaPower1 = mean(TFD1(:,ITheta),2);
    ThetaPower16 = mean(TFD16(:,ITheta),2);
    % Percentage power
    NormThetaPower1 = ThetaPower1/sum(ThetaPower1);
    NormThetaPower16 = ThetaPower16/sum(ThetaPower16);

    plot(T,NormThetaPower1)
    hold on
    plot(T,NormThetaPower16)
    xlabel('Time (s)')
    ylabel('\Theta Power')
    hold off

subplot(312)
    GamaPower1 = mean(TFD1(:,IGama),2);
    GamaPower16 = mean(TFD16(:,IGama),2);
    % Percentage power
    NormGamaPower1 = GamaPower1/sum(GamaPower1);
    NormGamaPower16 = GamaPower16/sum(GamaPower16);

    plot(T,NormGamaPower1)
    hold on
    plot(T,NormGamaPower16)
    xlabel('Time (s)')
    ylabel('Gama Power')
    hold off

subplot(313)
    HFOPower1 = mean(TFD1(:,IHFO),2);
    HFOPower16 = mean(TFD16(:,IHFO),2);
    % Percentage power
    NormHFOPower1 = HFOPower1/sum(HFOPower1);
    NormHFOPower16 = HFOPower16/sum(HFOPower16);

    plot(T,NormHFOPower1)
    hold on
    plot(T,NormHFOPower16)
    xlabel('Time (s)')
    ylabel('HFO Power')
    hold off
    
% %%%% SAME AS ABOVE BUT NOT NORMALIZED
% figure(9)
% subplot(311)
%     ThetaPower1 = mean(TFD1(:,ITheta),2);
%     ThetaPower16 = mean(TFD16(:,ITheta),2);
%     
%     plot(T,ThetaPower1)
%     hold on
%     plot(T,ThetaPower16)
%     xlabel('Time (s)')
%     ylabel('\Theta Power')
%     hold off
% 
% subplot(312)
%     GamaPower1 = mean(TFD1(:,IGama),2);
%     GamaPower16 = mean(TFD16(:,IGama),2);
%     
%     plot(T,GamaPower1)
%     hold on
%     plot(T,GamaPower16)
%     xlabel('Time (s)')
%     ylabel('Gama Power')
%     hold off
% 
% subplot(313)
%     HFOPower1 = mean(TFD1(:,IHFO),2);
%     HFOPower16 = mean(TFD16(:,IHFO),2);
%     
%     plot(T,HFOPower1)
%     hold on
%     plot(T,HFOPower16)
%     xlabel('Time (s)')
%     ylabel('HFO Power')
%     hold off


%% Question 10
% 10) Plot, in different subplots, scatter plots
% of all possible combinations of the time series
% in 9 (eg, theta ch 1 vs theta ch 16; theta ch 1
% vs gamma ch 16, etc), indicating as a title the
% linear correlation between the analyzed series.

ThetaPower1;
GamaPower1;
NormHFOPower1;

figure(10)
subplot(5,3,1)
    plot(ThetaPower1,GamaPower1,'ko')
    xlabel('\theta power electrode 1')
    ylabel('gama power electrode 1')
    title('Theta Ch 1 vs Gama Ch 1')

subplot(5,3,2)
    plot(ThetaPower1,HFOPower1,'ko')
    xlabel('\theta power electrode 1')
    ylabel('\HFO power electrode 1')
    title('Theta Ch 1 vs HFO Ch 1')
    
subplot(5,3,3)
    plot(ThetaPower1,ThetaPower16,'ko')
    xlabel('\theta power electrode 1')
    ylabel('\theta power electrode 16')
    title('Theta Ch 1 vs Theta Ch 16')
    
subplot(5,3,4)
    plot(ThetaPower1,GamaPower16,'ko')
    xlabel('\theta power electrode 1')
    ylabel('gama power electrode 16')
    title('Theta Ch 1 vs Gama Ch 16')
    
subplot(5,3,5)
    plot(ThetaPower1,HFOPower16,'ko')
    xlabel('\theta power electrode 1')
    ylabel('HFO power electrode 16')
    title('Theta Ch 1 vs HFO Ch 16')
    
    
subplot(5,3,6)
    plot(GamaPower1,HFOPower1,'ko')
    xlabel('gama power electrode 1')
    ylabel('HFO power electrode 1')
    title('Gama Ch 1 vs HFO Ch 1')
    
subplot(5,3,7)
    plot(GamaPower1,ThetaPower16,'ko')
    xlabel('gama power electrode 1')
    ylabel('\theta power electrode 16')
    title('Gama Ch 1 vs Theta Ch 16')
    
subplot(5,3,8)
    plot(GamaPower1,GamaPower16,'ko')
    xlabel('gama power electrode 1')
    ylabel('gama power electrode 16')
    title('Gama Ch 1 vs Gama Ch 16')

subplot(5,3,9)
    plot(GamaPower1,HFOPower16,'ko')
    xlabel('gama power electrode 1')
    ylabel('HFO power electrode 16')
    title('Gama Ch 1 vs HFO Ch 16')
    

subplot(5,3,10)
    plot(HFOPower1,ThetaPower16,'ko')
    xlabel('HFO power electrode 1')
    ylabel('\theta power electrode 16')
    title('HFO Ch 1 vs Theta Ch 16')
    
subplot(5,3,11)
    plot(HFOPower1,GamaPower16,'ko')
    xlabel('HFO power electrode 1')
    ylabel('gama power electrode 16')
    title('HFO Ch 1 vs Gama Ch 16')
    
subplot(5,3,12)
    plot(HFOPower1,HFOPower16,'ko')
    xlabel('HFO power electrode 1')
    ylabel('HFO power electrode 16')
    title('HFO Ch 1 vs HFO Ch 16')
    
    
subplot(5,3,13)
    plot(ThetaPower16,GamaPower16,'ko')
    xlabel('\theta power electrode 16')
    ylabel('gama power electrode 16')
    title('Theta Ch 16 vs Gama Ch 16')
    
subplot(5,3,14)
    plot(ThetaPower16,HFOPower16,'ko')
    xlabel('\theta power electrode 16')
    ylabel('HFO power electrode 16')
    title('Theta Ch 16 vs HFO Ch 16')
    
    
subplot(5,3,15)
    plot(GamaPower16,HFOPower16,'ko')
    xlabel('gama power electrode 16')
    ylabel('HFO power electrode 16')
    title('Gama Ch 16 vs HFO Ch 16')
set( findall(10, '-property', 'fontsize'), 'fontsize', 5)
    
    
%% Question 11
% 11) Compute and plot coherence spectra of
% channels 2 to 16 in relation to channel 1.  

clc, clf, clear, close all
load('LFPprobe.mat')

dt = 1/srate; % in s
Tmax = length(LFPprobe)/srate; % 600 s (10 min)
t = dt:dt:Tmax; % in s

windowlength = 1*srate;
overlap = [];
nfft = 2^16;

% possible pré alocation for speed. Gain speed of 7% 
CAll = zeros(16,32769); % Elapsed time is 83.6s.
% CAll = []; %Elapsed time is 89.4s.

tic
for i = 2:16
i;
[Cxy, F] = mscohere(LFPprobe(1,:),LFPprobe(i,:), ...
    windowlength,overlap,nfft,srate);
CAll(i,:) = Cxy;
end
toc

figure(11)
plot(F,CAll(2:16,:))
xlim([0 max(F)])
xlabel('Frequency (Hz)')
ylabel('Coherence')
legend('Ch1-Ch2','Ch1-Ch3','Ch1-Ch4',...
    'Ch1-Ch5','Ch1-Ch6','Ch1-Ch7', ...
    'Ch1-Ch8','Ch1-Ch9','Ch1-Ch10',...
    'Ch1-Ch11','Ch1-Ch12','Ch1-Ch13',... 
    'Ch1-Ch14','Ch1-Ch15','Ch1-Ch16', 'FontSize',6)
legend('Location','southeast')
legend('boxoff')
set( findall(11, '-property', 'fontsize'), 'fontsize', 5)


%% Question 12
% 12) Compute the average coherence in the theta
% band and also in the gamma band for each pair of
% channels analyzed in 11.  

ITheta = find(F(:,1) > 5 & F(:,1)<10);
IGama = find(F(:,1) > 30 & F(:,1)<100);

CAllTheta = mean(CAll(:,ITheta),2);
CAllGama = mean(CAll(:,IGama),2);


%% Question 13
% 13) Plot as a line graph the values obtained in
% 12, indicating on the Y axis the location
% (“height”) of the electrode (consider 100 to
% 1600 ?m) and on the X axis the average coherence
% value.    

altura = 100:100:1600;

figure (13)
plot(CAllTheta(2:end),altura(2:end),'k^-','MarkerFaceColor','k')
hold on
plot(CAllGama(2:end),altura(2:end),'r^-','MarkerFaceColor','r')
hold off
axis ij
xlabel('Cohrence')
ylabel('Profundidade (\mum)')
legend('\Theta Cohrence','Gama Cohrence')
legend('Location','northwest')
legend('boxoff')

%% Question 14
% 14) Compute and plot the coherogram between channels 1 and 16.

clc, clf, clear, close all
load('LFPprobe.mat')

dt = 1/srate;  % in s
Tmax = length(LFPprobe)/srate;  % in s
t = dt:dt:Tmax;  % in s

win = 6*srate; % window for Cxy computation
step = 0.1*win;
Nwin = (length(t)-win)/step+1;

cohwin = 2*srate;
nfft = 2^16;
overlap = [];

for nwin = 1:Nwin
nwin
winidx = (1:win) + (nwin-1)*step; 
[Cxy, F] = mscohere(LFPprobe(1,winidx),...
        LFPprobe(16,winidx),cohwin,overlap,nfft,srate);
Coherogram(nwin,:) = Cxy;
T(nwin) = t(winidx(win/2));    
end

clf
imagesc(T,F,Coherogram')
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

h = colorbar;
h.Label.String = 'Coherence';

%% Question 15
% 15) Compute and plot using imagesc the average
% coherence in a given frequency band (eg theta or
% gamma) for all possible combinations of channel
% pairs   

clc, clf, clear, close all
load('LFPprobe.mat')

dt = 1/srate;  % in s
Tmax = length(LFPprobe)/srate;  % in s
t = dt:dt:Tmax;  % in s

% Calculating the Coherence Spectrum
window = 5*srate; % window size in seconds
overlap = []; % window/2; %3*rate; % size of overlap between windows
nyquist = srate/2;

crosscoherenceTheta = zeros(16);
crosscoherenceGama = zeros(16);
for i = 1:16
    i
    SumTheta = zeros(16,1);
    SumGama = zeros(16,1);
    for j = i:16
        [Cxy, F] = mscohere(LFPprobe(i,:),LFPprobe(j,:),...
                            window,overlap,2^14,srate);
        ITheta = find(F(:,1) > 5 & F(:,1)<10);
        IGama = find(F(:,1) > 30 & F(:,1)<100);
        SumTheta(j) = sum(Cxy(ITheta));
        SumGama(j) = sum(Cxy(IGama));
    end
    crosscoherenceTheta(i,:) = SumTheta;
    crosscoherenceGama(i,:) = SumGama;
end

[macs,~] = max(abs(crosscoherenceTheta),[],2);
m = max(macs);
clear macs
ncrosscoherenceTheta = crosscoherenceTheta/m;

[macs,~] = max(abs(crosscoherenceGama),[],2);
m = max(macs);
clear macs
ncrosscoherenceGama = crosscoherenceGama/m;

figure(15)
subplot(1,2,1)
image(ncrosscoherenceTheta,'CDataMapping','scaled')
xlabel('Ch')
ylabel('Ch')
colorbar
h = colorbar;
h.Label.String = 'Theta Coherence';

subplot(1,2,2)
image(ncrosscoherenceGama,'CDataMapping','scaled')
xlabel('Ch')
ylabel('Ch')
colorbar
h = colorbar;
h.Label.String = 'Gama Coherence';











%%%%%%%%%%%%%%%% End %%%%%%%%%%%%%%%%%%%%%%