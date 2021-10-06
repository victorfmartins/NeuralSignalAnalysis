%% Lesson 3 = 09/03/2021

%% Fourier
% creates a composite sine signal and computes its
% Fourier coef at a differente freq

% setup
clc, clf, clear
srate = 1000; f = 50; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in s
X = sin(2*pi*f*t)+sin(2*pi*15*t); % compound sine signal

% plots the signal
subplot(3,1,1)
    plot(t,X)
    xlabel('Time (s)')
    ylabel('mV')
    title('LFP')

% plots 
subplot(3,1,[2 3])
    ff = 50; % freq Fourier Kernell
    K = exp(-1i*2*pi*ff*t); % Kernell
    FX = mean(X.*K); % Fourier coefficient on Kernell frequency (ff)
    Power = FX*conj(FX); % power at Kernell freq(ff)

    % plots lenght(t) points
    plot(X.*K,'yo')
    hold on
    
    %plots Fourier coef (vector)
    plot([0 real(FX)],[0 imag(FX)],'r-','linew',2)
    plot(real(FX),imag(FX),'ro','markerf','r','markers',7)
    hold off
    axis square
    xlim([-1 1])
    ylim(xlim())
    title(['Power at ' num2str(ff) ' Hz = ' num2str(Power)])

%% Power Spectrum
% computa o PSD sobre um vetor de frequencias
% o PSD é o poder do sinal across varias frequencias

% setup
srate = 1000; f = 15; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in s
X = sin(2*pi*f*t)+0*sin(2*pi*15*t); % sinal senoidal composto

% determina as freq para as quais obteremos seu poder.
freq = 0:0.01:80;
clear PSD

% computa e guarda o poder em cada frequencia - PSD
count = 0;
for ff = freq % freq Fourier Kernell
count = count+1;
K = exp(-1i*2*pi*ff*t);
FX = sum(X.*K);
PSD(count) = FX*conj(FX)/Tmax;
end

% plota o PSD
subplot(1,1,1)
    plot(freq,PSD,'k-'); hold on
    % cria um vetor booleano com apena um valor True e encontra sua posição
    % plota sua posição
    I = find(freq==15.1); 
    plot(freq(I),PSD(I),'ro','markerf','m'); hold off
    xlabel('Frequency (Hz)')
    ylabel('Power')
    xlim([10 20])
    title('Power Spectrum')

dF = freq(2)-freq(1);

%% Hamming window and Hann window
% Hamming: bell window that does not go to zero at the ends
% Hann: bell window that goes to zero at the ends

clf

% setup
srate = 1000; f = 15; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in s
X = sin(2*pi*f*t); % creates and plots a sine signal

subplot(3,1,1)
    plot(t,X)
    xlabel('Time (s)')
    ylabel('mV')
    title('LFP')

subplot(3,1,2)
    % bell window that goes to zero at the ends
    W1 = hann(length(X))'; 
    % W1 = rectwin(length(X))'; % unit window with zero at the ends
    plot(t,W1)
    % bell window that does not go to zero at the ends
    W2 = hamming(length(X))'; 
    hold on
    plot(t,W2)
    hold off
    legend('hann', 'hanning')
    legend('boxoff')
    xlabel('Time (s)')
    title('hann and hanning window')

subplot(3,1,3)
    plot(t,X.*W1)
    hold on
    plot(t,X.*W2)
    hold off
    legend('hann', 'hanning')
    legend('boxoff')
    xlabel('Time (s)')
    title('signal-hann and hanning window dot product')

%% PSDs with Hann/Hanning Win

freq = 0:0.01:80;
clear PSD PSDW1 PSDW2


% same as in PSD but multiplying by dt (why?) and windowing (W1 and W2)
% because this is the PDS formula, it is an integral
count = 0;
for ff = freq % freq Fourier Kernell
count = count+1;
K = exp(-1i*2*pi*ff*t);
FX = sum(X.*K)*dt;
PSD(count) = FX*conj(FX)/Tmax;

FX = sum((X.*W1).*K)*dt*1.63; %length(W1)/sum(W1);
PSDW1(count) = FX*conj(FX)/Tmax;

FX = sum((X.*W2).*K)*dt*1.59; %length(W2)/sum(W2);
PSDW2(count) = FX*conj(FX)/Tmax;

end

% plot PSDs
subplot(313)
    plot(freq,PSD,'k-')
    hold on
    plot(freq,PSDW1,'b-')
    plot(freq,PSDW2,'r-')
    % plot frequency marker 13.5
    I = find(freq==13.5);
    plot(freq(I),PSD(I),'ro','markerf','m')
    hold off
    xlabel('Frequency (Hz)')
    ylabel('Power')
    xlim([10 20])
    legend('PSD', 'PSD with Hann', 'PSD with Hanning')
    legend('boxoff')
    title('Comparison')


%% Usando a funcao built in do Matlab
% using pwelch(.)
% Pxx is the PSD
% pwelch uses Hanning Win

% setup
srate = 1000; f = 15; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in s
X = sin(2*pi*f*t);

window = length(X);
overlap = []; % defaut = 50%

% With rectwin(.) the side lopes appear (where do they come from??)
% % % % nfft = 2^16; % since srate = 1000 nfft default would be 2^10. 2^16
% just to make the graph more smoth
% % % % % using retangular window
% % % % [Pxx,F]=pwelch(X,rectwin(window),overlap,nfft,srate);
% % % % freqresolution = srate/nfft

% defining a frequency vector
freq = 0:0.01:80;

% freq = [0,1:5,6:0.5:14,14:0.01:20,25:5:80];
[Pxx,F]=pwelch(X,window,overlap,freq,srate);
% nfft default (padrao) e a menor potencia de 2 
% que seja maior que o tamanho da janela

subplot(313)
    hold on
    plot(F,Pxx,'y-','linew',3)
    xlabel('Frequency (Hz)')
    ylabel('Power (mV^2/Hz)')
    hold off
    % xlim([0 20])


%% the uncertainty principle in pwelch
% in stationary signals:
% the greater the length(win) the greater the
% confidence in the frequency and lower confidence
% in time 

% setup
srate = 1000; f = 5; % in Hz
dt = 1/srate; Tmax = 10; t = dt:dt:Tmax; % in s
X = sin(2*pi*f*t);

% window = length(X)
window = 4*srate; % range from 1 to 4 to see the effect of Unc. Princ.
overlap = []; % default: overlap  = window/2

% nfft is the numerical resolution of freq. how many frequency do i want = nfft/2
nfft = 10000; % (using power of 2 makes the algorithm run faster)
freqresolution = srate/nfft; % is the distance in freq that one freq has from the other

[Pxx,F]=pwelch(X,window,overlap,nfft,srate);

clf
subplot(111)
    plot(F,Pxx,'k-','linew',1)
    xlabel('Frequency (Hz)')
    ylabel('Power (mV^2/Hz)')
    xlim([0 20])
    title('Uncertainty Principle')








