%% Q1 a Q6 = Aula 09/03/2021; Q6 a Q11 = Aula 11/03/2021

clear, clc, clf
load('LFP_HG_HFO.mat')

%% Q1 a Q6 = Aula 09/03

srate = 1000; % in Hz (1/s)
dt = 1/srate; % in s
Tmax = length(lfpHG)*dt; % 300 s (5 min)
tvector = dt:dt:Tmax;
% tvector = (1:length(lfpHG))*dt;

plot(tvector,lfpHG)
hold on
plot(tvector,lfpHFO-1.2)
hold off

hold on
plot([11.5 11.5],[-1.5 -2]+0.1,'k-','linewidth',2)
plot([11.5 12],[-1.9 -1.9],'k-','linewidth',2)

text(11.6,-1.95,'0.5 s')

text(11.2,-1.8,'0.5 mV')


axis off

hold off

xlim([10 15])

xlabel('Time (s)')
ylabel('mV')


%% Q7

% % Compute a densidade espectral de pot?ncia (PSD) de cada canal 
% % usando o m?todo do Welch, com janelas de 2 segundos, 
% % 50% de sobreposi??o, e resolu??o num?rica de frequ?ncias 
% % de no m?nimo 0.1 Hz.

% detrend  <<<<<<<

%nfft = 10000;
nfft = 2^14;

% nfft = 16383; % esse é o teste que mostra que potencia de 2 vai mais ráp.
% nfft = 16384;

tic
[PxxHG,~] = pwelch(lfpHG,2*srate,1*srate,nfft,srate);
[PxxHFO,F] = pwelch(lfpHFO,2*srate,1*srate,nfft,srate);
toc

%% Q8

subplot(131)
plot(F,PxxHG)
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power mV^2/Hz')
title('LFP HG')

subplot(132)
plot(F,PxxHFO)
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power mV^2/Hz')
title('LFP HFO')

subplot(133)
plot(F,PxxHG)
hold on
plot(F,PxxHFO)
hold off
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power mV^2/Hz')
legend('lfpHG','lfpHFO')
legend('boxoff')

%% Q9 e Q10

I = find(F > 5 & F < 10);

ThetaMHG = mean(PxxHG(I));
ThetaMHFO = mean(PxxHFO(I));

% % subplot(132)
% % plot(F,PxxHFO)
% % hold on
% % plot(F(I),PxxHFO(I),'y-')
% % hold off
% % 
% % xlim([0 20])
% % xlabel('Frequency (Hz)')
% % ylabel('Power mV^2/Hz')
% % title('LFP HFO')

subplot(131)

bar([1 2],[ThetaMHG ThetaMHFO],0.5,'c')
ylabel('Theta power')

label{1} = 'lfpHG';
label{2} = 'lfpHFO';

set(gca,'xtick',[1 2],'xticklabel',label)
set(gcf,'color','w') %[0.5 0.5 0.5])

%% Q11

subplot(211)
plot(F,PxxHG)
hold on
plot(F,PxxHFO)
hold off
xlim([0 200])
xlabel('Frequency (Hz)')
ylabel('Power mV^2/Hz')

subplot(212)
plot(F,10*log10(PxxHG))
hold on
plot(F,10*log10(PxxHFO))
hold off
xlim([0 200])

xlabel('Frequency (Hz)')
% é p 10*log10(.) e não tem unidade pois é um ratio entre dois objetos 
% de mesmas dimensões. no nosso caso é entre Pxx/(1*mV^2/Hz)
ylabel('Power dB') 

%%


f = 0.001:0.001:0.1;

Pxx = 1./x;
plot(f,Pxx)


%%


x  = 0.001:0.001:100;

plot(x,log10(x))
hold on
plot([0 100],[0 0])
hold off
















