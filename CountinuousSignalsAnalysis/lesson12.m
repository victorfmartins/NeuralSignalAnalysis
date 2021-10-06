%% Lesson 12 = 15/04/2021; Wavelets
% building a Complex Morlet Wavelet 

clf, clc, clear
srate = 1000; ff = 15; % in Hz
dt = 1/srate; t_psi = -2:dt:2;

senoide1 = cos(2*pi*ff*t_psi);
senoide2 = sin(2*pi*ff*t_psi);

% parameters of a Gaussian function 
a = 1; % higth of the pico
m = 0; % position central (mean)
sd = 0.1; % standard deviation ("width")

% sd - on practice - determines the number of
% cycles has a formula relating these parameters  
NumeroCiclos = 7;
sd = NumeroCiclos/(2*pi*ff);

Gaussiana = a*exp(-(t_psi-m).^2/(2*sd^2));

PsiR = Gaussiana.*senoide1;
PsiI = Gaussiana.*senoide2;
Psi = PsiR + PsiI*1i;
PsiC = Gaussiana.*exp(1i*2*pi*ff*t_psi);

subplot(311)
plot(t_psi,senoide1,'b-')
hold on
plot(t_psi,senoide2,'r-')
plot(t_psi,Gaussiana,'k-')
hold off
xlabel('Time (s)')
xlim([-1 1])

subplot(312)
plot(t_psi,PsiR,'b-')
hold on
plot(t_psi,PsiI,'r-')
% plot(t_psi,real(PsiC),'k-')
% plot(t_psi,imag(PsiC),'k-')
hold off
xlabel('Time (s)')
xlim([-1 1])

subplot(313)
plot3(t_psi,PsiR,PsiI)
xlabel('Time (s)')
ylabel('Real part')
zlabel('Imaginary part')
xlim([-1 1])

%% Convoluting a signal with Complex Morlet Wavelet 

srate = 1000; f = 10; ff = 10; % in Hz
dt = 1/srate; t_psi = -0.3:dt:0.3; t_lfp = 0:dt:2;

NumeroCiclos = 5;
sd = NumeroCiclos/(2*pi*ff);
Psi = exp(-t_psi.^2/(2*sd^2)).*exp(1i*2*pi*ff*t_psi);
PsiR = real(Psi);
PsiI = imag(Psi);

LFP = sin(2*pi*f*t_lfp) + 0.5*sin(2*pi*(f+15)*t_lfp);
LFP = sin(2*pi*f*t_lfp);
LFP(1:1000) = 0;
LFP = LFP+0.5*sin(2*pi*(f+15)*t_lfp);

% LFP = [zeros(size(Psi)) LFP zeros(size(Psi))];
% t_lfp = dt*(1:length(LFP));

for j = 0:1:length(LFP)-length(Psi)

subplot(211)
plot(t_lfp,LFP)
hold on
plot(t_lfp((1:length(Psi))+j),PsiR,'g','linew',2)
plot(t_lfp((1:length(Psi))+j),PsiI,'r','linew',2) 
hold off

ylim([-1.5 1.5])
xlim([0 t_lfp(end)])

ConvR = sum(flip(PsiR).*LFP([1:length(Psi)]+j));
ConvR = ConvR/sum(PsiR.^2);

ConvI = sum(flip(PsiI).*LFP([1:length(Psi)]+j));
ConvI = ConvI/sum(PsiI.^2);

subplot(212)
plot(t_lfp(floor(length(Psi)/2)+j),ConvR,'g.')
hold on
plot(t_lfp(floor(length(Psi)/2)+j),ConvI,'r.')

ylim([-1.5 1.5])
xlim([0 t_lfp(end)])
pause(0.0001)
end

hold off

% % using the build in convolutioni

Conv = conv(LFP,Psi/(sum(PsiR.^2)),'same');
AmpEnv = abs(Conv);

subplot(212)
hold on
plot(t_lfp,real(Conv),'g-')
plot(t_lfp,imag(Conv),'r-')
plot(t_lfp,AmpEnv,'k-','linew',3)
hold off

%% Using Matlab's built in functions for Wavelets 

waveinfo('cmor')

% a complex morlet has two parameters:

% Fc - frequency of the senoid
% 
% Fb - bandwidth parameter

% Note that Fb = 2*sd^2

% defining Fb based on the number of cycles
NumeroCiclos = 5;
sd = NumeroCiclos/(2*pi*ff);
Fb = 2*sd^2;

%%

% sintaxe for construct a wavelet:

% wname = 'cmorFb-Fc'
wname = 'cmor0.5-2.3'

[MotherPsi,ttt] = wavefun(wname,15);
plot(ttt,real(MotherPsi),'g-')
hold on
plot(ttt,imag(MotherPsi),'r-')
plot([1 1],[-0.5 0.5])
hold off

% geting the center frequancy
F = centfrq(wname)

% wavefun(wname,15,0) % plot the caracteristic wavelet

%%

Fc = 10;
NumeroCiclos = 5;

sd = NumeroCiclos/(2*pi*Fc);
Fb = 2*sd^2;
wname = ['cmor' num2str(Fb) '-' num2str(Fc)]


[MotherPsi,ttt] = wavefun(wname,15);


plot(ttt,real(MotherPsi),'g-','linew',3)
hold on
plot(ttt,imag(MotherPsi),'r-','linew',3)
hold off

% xlim([-0.5 0.5])

%% Doing a continuos wavelet transform
% using the built in function of matlab

t_lfp = 0:dt:2;

% LFP = sin(2*pi*10*t_lfp) + 0.5*sin(2*pi*25*t_lfp);

% LFP = sin(2*pi*10*t_lfp);
% LFP(1:1000) = 0;
% LFP = LFP+0.5*sin(2*pi*25*t_lfp);

LFP = sin(2*pi*10*t_lfp)+sin(2*pi*40*t_lfp);
LFP(1:1000)=0;

% WT stands for Wavelet Transform

% the scale factor determines the frequencies that 
% will be analyzed, being inversely proportional 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Freq = 40;
scale = Fc/(Freq*dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WT = cwt(LFP,scale,wname);
Amp = abs(WT);

subplot(211)
plot(t_lfp,LFP)

subplot(212)
plot(t_lfp,real(WT),'g')
hold on
plot(t_lfp,imag(WT),'r')
plot(t_lfp,Amp,'k-','linew',2)
hold off

ylim([-30 30])

%% TFD continuos using Wavelets
srate = 1000;
dt = 1/srate;
Fc = 1;
NumeroCiclos = 7;


sd = NumeroCiclos/(2*pi*Fc);
Fb = 2*sd^2;
wname = ['cmor' num2str(Fb) '-' num2str(Fc)];
Freq = 1:0.5:100;
% scale = centfreq(wname)./(Freq*dt);
scale = Fc./(Freq*dt);

tic
WT = cwt(LFP,scale,wname);
toc

TFD = abs(WT);

clf

imagesc(t_lfp,Freq,TFD)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% Compraring with fourier

subplot(211)
imagesc(t_lfp,Freq,TFD)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar


subplot(212)
tic
[S,F,T,TFD2]=spectrogram(LFP,0.2*srate,[],Freq,srate);
toc
imagesc(T,F,TFD2)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar
















