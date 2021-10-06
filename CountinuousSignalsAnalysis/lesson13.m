%% Lesson 13 = 20/04/2021; Fase e Freq Estantania, PLV
%% Instantaneous Phase

% Setup
clear, clc, clf
srate = 1000; f = 7;
dt = 1/srate; t = dt:dt:2;

% Signal of freq f modulated by A (a sine of freq 1Hz)
A = 2+sin(2*pi*1*t);
LFP = A.*sin(2*pi*f*t);

%signal to complex signal to Envelope of Amplitude
LFPanalitico = hilbert(LFP);
HT = imag(LFPanalitico);
EnvAmp = abs(LFPanalitico);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phase = angle(LFPanalitico); % given in radians from -pi to pi
% Phase = Phase+pi;
% Phase(Phase<0) = Phase(Phase<0)+2*pi;
Phase = rad2deg(Phase);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:4:2000
h1 = subplot(311);
    plot(t,LFP,'b-')
    hold on
    plot(t,HT,'r-')
    plot(t,EnvAmp,'k','linew',2)
    plot(t(j),LFP(j),'bo','markerf','b')
    plot(t(j),HT(j),'ro','markerf','r')
    plot(t(j),EnvAmp(j),'ko','markerf','k')
    hold off
    xlabel('Time (s)')
    ylabel('mV')

h2 = subplot(312);
    plot(t,Phase,'g.')
    hold on
    plot(t(j),Phase(j),'go','markerf','g')
    hold off
    % ylim([0 2*pi])
    % ylim([0 360])
    % set(gca,'ytick',0:90:360)
    set(gca,'ytick',-180:90:180)
    title(['Inst Phase = ' num2str(Phase(j)) '^o'])
    xlabel('Time (s)')
    ylabel('Phase (^o)')

subplot(313)
    plot([0 LFP(j)],[0 HT(j)],'k-')
    hold on
    plot(LFP(j),HT(j),'ko','markerf','k')
    plot(LFP(j),0,'bo','markerf','b')
    plot(0,HT(j),'ro','markerf','r')
    hold off
    xlim([-3 3])
    ylim([-3 3])
    axis square

pause(0.001)
end

linkaxes([h1 h2],'x')

%% Instantaneous Frequency

% Setup
clear, clc, clf
srate = 1000; f = 7;
dt = 1/srate; t = dt:dt:4;

% Modulation (changing) of Amplitude 
A = 2+sin(2*pi*1*t);
LFP = A.*sin(2*pi*f*t);

% Frequency modulation done wrong
% F = 2*pi*10*t;
% freq = 10+4*sin(2*pi*1*t);
% freq = 10 + randn(size(t));
% LFP = A.*sin(2*pi*freq.*t);

% Modulation (changing) of Frequency and Amplitude 
% programming frequency modulation properly
freq = 10+4*sin(2*pi*1*t);
LFP = A.*sin(cumsum(2*pi*freq.*dt));

LFPanalitico = hilbert(LFP);
HT = imag(LFPanalitico);
EnvAmp = abs(LFPanalitico);

Phase = angle(hilbert(LFP));
% % % Freq = diff(Phase)/(2*pi*dt); % wrong formula
PhaseUnwrapped = unwrap(Phase);
Freq = diff(PhaseUnwrapped)/(2*pi*dt); 
% alternative form
Freq2 = angle(exp(1i*diff(Phase)))/(2*pi*dt);


h1 = subplot(311);
    plot(t,LFP,'b-')
    hold on
    plot(t,HT,'r-')
    plot(t,EnvAmp,'k','linew',2)
    hold off
    xlabel('Time (s)')
    ylabel('mV')

h2 = subplot(312);
    plot(t,Phase,'g.')
    % hold on
    % plot(t,PhaseUnwrapped,'g-')
    % hold off
    xlabel('Time (s)')
    ylabel('Phase (rad)')

h3 = subplot(313);
    plot(t(1:end-1)+dt/2,Freq)
    hold on
    plot(t(1:end-1)+dt/2,Freq2)
    hold off
    xlabel('Time (s)')
    ylabel('Inst Freq (Hz)')

linkaxes([h1 h2 h3],'x')

%% Phase-locking metrics

% Setup
clear, clc, clf
srate = 1000; f = 8;
dt = 1/srate; t = dt:dt:4;

ruido1 = 3*randn(size(t));
ruido2 = 3*randn(size(t));

LFP1 = sin(2*pi*f*t) + ruido1;
LFP2 = sin(2*pi*f*t+pi/2) + ruido2;

LFP1filtrado = eegfilt(LFP1,srate,6,10);
LFP2filtrado = eegfilt(LFP2,srate,6,10);

% LFP1filtrado = eegfilt(LFP1,srate,50,60);
% LFP2filtrado = eegfilt(LFP2,srate,50,60);

Phase1 = angle(hilbert(LFP1filtrado));
Phase2 = angle(hilbert(LFP2filtrado));

% DeltaPhase = Phase2-Phase1; % > can yeild problems
DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % recomended method

%%% 
figure(1)
h1 = subplot(311);
    plot(t,LFP1,'b')
    hold on
    plot(t,LFP2,'r')
    plot(t,LFP1filtrado-3,'b')
    plot(t,LFP2filtrado-3,'r')
    hold off
    xlim([2 3])

h2 = subplot(312);
    plot(t,Phase1,'b.')
    hold on
    plot(t,Phase2,'r.')
    hold off
    ylabel('\Phi (rad)')
    set(gca,'ytick',-pi:pi/2:pi,...
        'yticklabel',{'-pi','-pi/2','0','pi/2','pi'})
    xlim([2 3])

h3 = subplot(313);
    plot(t,rad2deg(DeltaPhase),'k.')
    ylabel('\Delta\Phi (^o)')
    xlabel('Time (s)')
    xlim([2 3])
    ylim([-180 180])
linkaxes([h1 h2 h3],'x')

%%% seeing phase difference histograms 
figure(2)
subplot(221)
    phasebins = -170:20:170;
    [counts, phasebins]=hist(DeltaPhase,deg2rad(phasebins));
    bar(rad2deg(phasebins),counts,'k')
    xlabel('\Delta\Phi (^o)')
    ylabel('counts')
    set(gca,'xtick',-180:90:180)

subplot(222)
    rose(DeltaPhase,18) % number of phase bins
    title('\Delta\Phi (^o)')

subplot(223)
    compass(exp(1i*DeltaPhase(1:20:end)))

subplot(224)
    polar(angle(exp(1i*DeltaPhase(1:20:end))),...
        abs(exp(1i*DeltaPhase(1:20:end))),'ko')
    VetorMedio = mean(exp(1i*DeltaPhase));
    hold on
    compass(VetorMedio,'r')
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computing the Phase-Locking Value
    PLV = abs(mean(exp(1i*DeltaPhase)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    title(['PLV = ' num2str(PLV)])











