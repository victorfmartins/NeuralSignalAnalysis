%% Lesson 14 = 22/04/2021; Especter of PLV, EPL, Kappa and ITC
% Especter of PLV

% Setup
clear, clc, clf
srate = 1000; f = 8;
dt = 1/srate; t = dt:dt:4;

ruido1 = 3*randn(size(t));
ruido2 = 3*randn(size(t));

LFP1 = sin(2*pi*f*t) + ruido1;
LFP2 = sin(2*pi*f*t+0*pi/2) + ruido2;

freqvector = 4:0.5:100; % in Hz
bandwidth = 4; % in Hz
count = 0;
PLVspectrum = zeros(size(freqvector));

h = waitbar(0,'Computing PLV Spectrum');
for f = freqvector
    count = count+1;
    waitbar(f/freqvector(end),h)

    LFP1filtrado = eegfilt(LFP1,srate,f,f+bandwidth);
    LFP2filtrado = eegfilt(LFP2,srate,f,f+bandwidth);    

    Phase1 = angle(hilbert(LFP1filtrado));
    Phase2 = angle(hilbert(LFP2filtrado));

    DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % metodo recomendado

    PLV = abs(mean(exp(1i*DeltaPhase)));
    PLVspectrum(count) = PLV;
end
close(h)

plot(freqvector+bandwidth/2,PLVspectrum)
hold on
plot(freqvector+bandwidth/2,PLVspectrum.^2)
hold off
xlabel('Frequency (Hz)')
ylabel('PLV')
legend('PLV','PLV^2')
title('PLV Spectrum')

%% Finite sample bias


clear, clc, clf
srate = 1000;
dt = 1/srate;
t = dt:dt:4;

ruido1 = 1*randn(size(t));
ruido2 = 1*randn(size(t));

LFP1 = sin(2*pi*8*t) + ruido1;
LFP2 = sin(2*pi*8*t+0*pi/2) + ruido2;

LFP1filtrado = eegfilt(LFP1,srate,60,70);
LFP2filtrado = eegfilt(LFP2,srate,60,70);

LengthVector = 5:5:4000;
count = 0;
PLV_epochlength = size(LengthVector);

for EpochLength = LengthVector; % in seconds
    count = count+1;

    Phase1 = angle(hilbert(LFP1filtrado(1:EpochLength)));
    Phase2 = angle(hilbert(LFP2filtrado(1:EpochLength)));

    DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % recomended method
    PLV = abs(mean(exp(1i*DeltaPhase)));
    PLV_epochlength(count) = PLV;
end

plot(LengthVector,PLV_epochlength)

%% other phase locking metrics

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:4;

ruido1 = 5*randn(size(t));
ruido2 = 5*randn(size(t));

LFP1 = sin(2*pi*8*t) + ruido1;
LFP2 = sin(2*pi*8*t+0*pi/2) + ruido2;

LFP1filtrado = eegfilt(LFP1,srate,60,100);
LFP2filtrado = eegfilt(LFP2,srate,60,100);

Phase1 = angle(hilbert(LFP1filtrado));
Phase2 = angle(hilbert(LFP2filtrado));

DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % recomended method


% % seeing phase difference histograms 
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
% z1 = 1 + 1*1i;
% z2 = -1 + 1*1i;
% compass([z1 z2]);

compass(exp(1i*DeltaPhase(1:20:end)))


subplot(224)

polar(angle(exp(1i*DeltaPhase(1:20:end))),...
    abs(exp(1i*DeltaPhase(1:20:end))),'ko')

VetorMedio = mean(exp(1i*DeltaPhase));

hold on
compass(VetorMedio,'r')
hold off


% computing the Phase-Locking Value
PLV = abs(mean(exp(1i*DeltaPhase)));
title(['PLV = ' num2str(PLV)])

% % ENTROPY BASED INDEX
% Defining an entropy-based coupling index

% turning histogram counts into probabilities 

phasebins = -170:20:170;
[counts, phasebins]=hist(DeltaPhase,deg2rad(phasebins));

p = counts/sum(counts);
Nbins = length(p);

H = -sum(p(p>0).*log(p(p>0)));
Hmax = log(Nbins);

% entropy-based phase locking
EPL = (Hmax-H)/Hmax;

subplot(223)
bar(rad2deg(phasebins),p,'k')
xlabel('\Delta\Phi (^o)')
ylabel('probability')
set(gca,'xtick',-180:90:180)
title(['Entropy based index = ' num2str(EPL)])


% %% VON MISES DISTRIBUTION tutorial
% clf
% 
% % aka Distribuicao Gaussiana no Circulo
% meanphi = deg2rad(-20) % mean angle
% kappa = 100 % "inverse" of the standard deviation
% 
% phi = -pi:0.01:pi;
% 
% VonMises = exp(kappa*cos(phi-meanphi))/(2*pi*besseli(0,kappa));
% 
% plot(rad2deg(phi),VonMises,'y-','linewidth',3)
% xlabel('\phi (deg)')
% ylabel('Prob')
% % ylim([0 1])
% 
% % %


phasebins = -170:20:170;
[counts, phasebins]=hist(DeltaPhase,deg2rad(phasebins));


% from the Circular Statistics Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[meantheta,anglestrength,sigma,confangle,kappa]=anglemean(DeltaPhase); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = -pi:0.01:pi;
VonMises = exp(kappa*cos(phi-meantheta))/(2*pi*besseli(0,kappa));


p = counts/sum(counts);
Nbins = length(p);

dPhi = phasebins(2)-phasebins(1);

subplot(222)
bar(rad2deg(phasebins),p,'k')
hold on
plot(rad2deg(phi),VonMises*dPhi,'y-','linewidth',3)
hold off
xlabel('\Delta\Phi (^o)')
ylabel('probability')
set(gca,'xtick',-180:90:180)
title(['Kappa = ' num2str(kappa)])






%% Especter of Phase locking again

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:20;

ruido1 = 3*randn(size(t));
ruido2 = 3*randn(size(t));

LFP1 = sin(2*pi*8*t) + ruido1;
LFP2 = sin(2*pi*8*t+0*pi/2) + ruido2;

freqvector = 0:0.5:100; % in Hz
bandwidth = 4; % in Hz

clear PLVspectrum EPLspectrum Kappaspectrum

h = waitbar(0,'Computing PLV Spectrum');

count = 0;
for f = freqvector
    count = count+1;
waitbar(f/freqvector(end),h)

LFP1filtrado = eegfilt(LFP1,srate,f,f+bandwidth);
LFP2filtrado = eegfilt(LFP2,srate,f,f+bandwidth);    
    
Phase1 = angle(hilbert(LFP1filtrado));
Phase2 = angle(hilbert(LFP2filtrado));

DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % recomended method

PLV = abs(mean(exp(1i*DeltaPhase)));
PLVspectrum(count) = PLV;

phasebins = -170:20:170;
[counts, phasebins]=hist(DeltaPhase,deg2rad(phasebins));
p = counts/sum(counts);
Nbins = length(p);
H = -sum(p(p>0).*log(p(p>0)));
Hmax = log(Nbins);
% entropy-based phase locking
EPL = (Hmax-H)/Hmax;
EPLspectrum(count) = EPL;

% Von mises kappa
[meantheta,anglestrength,sigma,confangle,kappa]=anglemean(DeltaPhase); 
Kappaspectrum(count) = kappa;

end

close(h)

%%

subplot(221)
plot(freqvector+bandwidth/2,PLVspectrum)
hold on
plot(freqvector+bandwidth/2,PLVspectrum.^2)
hold off
xlabel('Frequency (Hz)')
ylabel('PLV')
legend('PLV','PLV^2')

subplot(222)
plot(freqvector+bandwidth/2,EPLspectrum)
xlabel('Frequency (Hz)')
ylabel('Entropy based index')


subplot(223)
plot(freqvector+bandwidth/2,Kappaspectrum)
xlabel('Frequency (Hz)')
ylabel('Kappa')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OSCILLATORY RESET

clf, clear, clc

figure(1)

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

for ntrial = 1:30
   
    LFP = 0.5*sin(2*pi*6*t + 2*pi*rand);
    LFP(1001:1200) = 0.5*sin(2*pi*6*t(1001:1200)+0.5*randn);
    
    Phase = angle(hilbert(LFP));
    PhaseAll(ntrial,:)=Phase;
    
    h1 = subplot(5,1,[1 4]);
%     subplot(211)
    plot(t-1,LFP+ntrial,'k-')
%     hold on
%     subplot(212)
%     plot(t-1,0.1*Phase+ntrial,'k.')
%     
    hold on
    
end

plot([0 0],[0 ntrial+1],'r--','linew',2)
hold off
ylabel('Trial #')
% xlabel('Time (s)')
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ITC = abs(mean(exp(1i*PhaseAll)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h2 = subplot(5,1,5);
plot(t-1,ITC,'k','linew',2)
xlabel('Time (s)')
ylabel('ITC')
ylim([0 1])

linkaxes([h1 h2],'x')

%%

for j = 1:10:2000 % time index
figure(2)
compass(exp(1i*PhaseAll(:,j)))
meanvector = mean(exp(1i*PhaseAll(:,j)));
hold on
compass(meanvector,'k-')
hold off
% ITC stands for Inter Trial Coherence
ITC = abs(meanvector);
title(['ITC at ' num2str(t(j)) ' s = ' num2str(ITC)])
pause(0.01)
end


%% We call it ITC (inter trial coherence) if we analyze
% the same channel on different trials.
% if we analyze different channels in the same
% trial, this same metric corresponds to a
% "global" synchrony index. GSI (global synchrony
% index) or GPC (global phase coherence)    

% ITC spectrum

clf, clear, clc

figure(1)

srate = 1000;
dt = 1/srate;
t = dt:dt:5;

freqvector = 0:2:60;

for ntrial = 1:100
    ntrial;
    LFP = 0.5*sin(2*pi*30*t + 2*pi*rand);
    LFP(2500:3000) = 0.5*sin(2*pi*30*t(2500:3000)+0.5*randn);
    LFP = LFP + randn(size(LFP));
    
    count = 0;
    for f = freqvector
    count = count+1;
    
    LFPfiltrado = eegfilt(LFP,srate,f,f+4);
    Phase = angle(hilbert(LFPfiltrado));
    PhaseAll(ntrial,count,:) = Phase;
   
    end
    
end

%%

ITCspectrum = abs(mean(exp(1i*PhaseAll(:,:,:))));
ITCspectrum = squeeze(ITCspectrum);

imagesc(t-2.5,freqvector + 2,ITCspectrum)
axis xy
ylabel('Frequency (Hz)')
xlabel('Peristimulus time (s)')
title('Inter trial coherence')
colorbar








