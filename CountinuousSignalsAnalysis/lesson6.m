%% Lesson 6 = 18/03/2021; (Lesson 5 generated list 2 resolution file)

%% Understanding Phase Coherence

%setup
clear, clf, clc,

% create complex x
R = 3;
theta = deg2rad(45);
Fx = R*exp(1i*theta);

% create complex y
R = 0.5;
theta = deg2rad(145);
Fy = R*exp(1i*theta);

% make normalization of x and y
normFx = Fx/abs(Fx);
normFy = Fy/abs(Fy);

subplot(111)
    % plot x and norm of x
    plot([real(Fx)],[imag(Fx)],'bo','markerf','b')
    hold on
    plot([0 real(Fx)],[0 imag(Fx)],'b-')
    plot([real(normFx)],[imag(normFx)],'bo','markerf','b')
    plot([0 real(normFx)],[0 imag(normFx)],'b-')

    % plot y and norm of y
    plot([real(Fy)],[imag(Fy)],'ro','markerf','r'); hold on
    plot([0 real(Fy)],[0 imag(Fy)],'r-')
    plot([real(normFy)],[imag(normFy)],'ro','markerf','r')
    plot([0 real(normFy)],[0 imag(normFy)],'r-')

    % plot unit circle
    plot(exp(1i*(0:0.01:2*pi)),'k-'); hold off

    % creates legend
    text([real(Fx)+0.1],[imag(Fx)],'Fx')
    text([real(Fy)+0.1],[imag(Fy)],'Fy')

    text([real(normFx)+0.1],[imag(normFx)],'Fx/|Fx|')
    text([real(normFy)+0.1],[imag(normFy)],'Fy/|Fy|')

    xlim([-3 3])
    ylim([-3 3])
    axis square
    title('vectors and their normalizations to the unitary circle')

%% Vector average of unit vectors as coherence preview
% creates N unit vectors and computes
% coherence is defined as the length of the mean vector of N vector

clear Vall
N = 20; % number of vectors

% plot N vectors
for nvetor = 1:N
% random angle unit vector
theta = pi/2 + 1*randn;
v = exp(1i*theta);

Vall(nvetor) = v;

plot([real(v)],[imag(v)],'bo','markerf','b')
hold on
plot([0 real(v)],[0 imag(v)],'b-')
plot(exp(1i*(0:0.01:2*pi)),'k-')

xlim([-1 1])
ylim([-1 1])
axis square
end

MeanVector = mean(Vall);
plot([real(MeanVector)],[imag(MeanVector)],'ko','markerf','k')
plot([0 real(MeanVector)],[0 imag(MeanVector)],'k-')
hold off

% coherence is the defined as the length of the
% mean vector (mean over all unitary vectors)
coherence = abs(MeanVector);

title(['Coherence = ' num2str(coherence)])


%% Coherence from the beginning
% generates two sinusoidal signals with the same
% frequency and a phase diff computes a conv at
% frequency ff for the entire signal (one window)
% of the conv of X and Y signal with K obtains the
% coherence point at frequency ff

%setup
clear, clf, clc,
srate = 1000;  f = 8; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in s

% create two signals of frequency f, with phase
% difference of phi and with noise 
phi = -deg2rad(45);
X = 2*sin(2*pi*f*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*f*t+phi)+0.3*randn(size(t));
    
ff = 8; % Frequency of Fourier Kernel
K = exp(-1i*2*pi*ff*t); % Kernel

h1 = subplot(311);
    plot(t,X,'b'), hold on
    plot(t,Y,'r')
    plot(t,real(K),'k'), hold off
    xlim([0 1])
    xlabel('Time (s)')
    ylabel('mv')
    title('Fake LFP')

XX = X.*K;
YY = Y.*K;
    
% actually, FX is the sum, we are using the 
% mean just for visualization
FX = mean(XX); disp(['FX = ' num2str(FX)])
FY = mean(YY); disp(['FY = ' num2str(FY)])

nFX = FX/abs(FX); disp(['nFX = ' num2str(nFX)])
nFY = FY/abs(FY); disp(['nFY = ' num2str(nFY)])

% plot to see how X.*K and Y.*K look like
h2 = subplot(312);
    plot(t,real(XX),'b','linew',2), hold on
    plot(t,real(YY),'r','linew',2)
%     plot(t,imag(XX),'b--')
%     plot(t,imag(YY),'r--')
    plot(t(end/2),real(nFX),'ko','MarkerFaceColor','b')
    plot(t(end/2),real(nFY),'ko','MarkerFaceColor','r')
    hold off
    xlim([0 1])
    xlabel('Time (s)')
    ylabel('mv')
    title('Pointwise Multiplication of Fake LFP with Kernel')

linkaxes([h1 h2])
% xlim([0 .2])
    
nFXY = nFX*conj(nFY); % coherence point of the freq ff
% nFXY = FX*conj(FY)/(abs(FX)*abs(FY));

subplot(3,1,3)
    plot([real(FX)],[imag(FX)],'bo')
    hold on
    plot([0 real(FX)],[0 imag(FX)],'b-')
    plot([real(FY)],[imag(FY)],'ro')
    plot([0 real(FY)],[0 imag(FY)],'r-')

    plot([real(nFX)],[imag(nFX)],'bo','markerf','b')
    plot([0 real(nFX)],[0 imag(nFX)],'b-')
    plot([real(nFY)],[imag(nFY)],'ro','markerf','r')
    plot([0 real(nFY)],[0 imag(nFY)],'r-')

    plot([real(nFXY)],[imag(nFXY)],'ko','markerf','k')
    plot([0 real(nFXY)],[0 imag(nFXY)],'k-')

    plot(exp(1i*(0:0.01:2*pi)),'k-')
    hold off

    xlim([-2 2])
    ylim([-2 2])
    axis square
    title('vectors and their normalizations to the unitary circle')

%% loop across (fake) time windows
% simulates N win each with two simple sine
% signals one with phase phi and the other not
% phi is constant in and over win
% computes the coherence point on freq ff for each win
% computes the coh over windows point on freq ff

%%% instead of calculating coh on one signal
%%% windows you can calculate on different small
%%% trials = intertrialcoh  

%%% running with ff = 7.5 and Tmax = 3 gives more
%%% coh than running with Tmax = 2. It was to be
%%% expected that with longer times (time of
%%% windows) the coh would reduce (which is vdd in
%%% general), but with Tmax = 3 it has beat and
%%% with Tmax = 2 it does not.
%%% Think about this when running your coherences


%setup
clear, clf, clc,
srate = 1000;  f = 8; ff = 7.5; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in s

% Fourier Kernel
K = exp(-1i*2*pi*ff*t);

clear nFXYAll
N = 100; % number of time windows

% compute and plot N. each point is the
% window-to-window coherence. each window has Tmax
% seconds and has a similar phase difference
for nwindow = 1:N

% create two signals with noise of frequency f,
% with phase difference phi and with noise 
phi = -deg2rad(0)+0.3*randn;
X = 3*sin(2*pi*8*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*8*t+phi)+0.3*randn(size(t));

% actually, FX is the sum, we are using the 
% mean just for visualization
FX = mean(X.*K);
FY = mean(Y.*K);

nFX = FX/abs(FX);
nFY = FY/abs(FY);

% coherence point of the freq $ff and win nwindow
nFXY = nFX*conj(nFY);

nFXYAll(nwindow) = nFXY;
% nFXY = FX*conj(FY)/(abs(FX)*abs(FY));

% These two blocks contain the instruction to
% print the nFx and nFy of each win 
% % plot([real(FX)],[imag(FX)],'bo')
% % hold on
% % plot([0 real(FX)],[0 imag(FX)],'b-')
% % plot([real(FY)],[imag(FY)],'ro')
% % plot([0 real(FY)],[0 imag(FY)],'r-')
% % 
% % plot([real(nFX)],[imag(nFX)],'bo','markerf','b')
% % plot([0 real(nFX)],[0 imag(nFX)],'b-')
% % plot([real(nFY)],[imag(nFY)],'ro','markerf','r')
% % plot([0 real(nFY)],[0 imag(nFY)],'r-')


plot([real(nFXY)],[imag(nFXY)],'ko','markerf','k')
hold on
plot([0 real(nFXY)],[0 imag(nFXY)],'k-')

plot(exp(1i*(0:0.01:2*pi)),'k-')
hold on

xlim([-2 2])
ylim([-2 2])
axis square

title(['Nwindow = ' num2str(nwindow)])
pause(0.005)
end

% nFXYAll gives us N coh points across N windows
% of Tmax seconds and frequency ff 
Cxy = abs(mean(nFXYAll));

text(0,1.5,['Coherence = ' num2str(Cxy)],'fontsize',16)

hold on

plot([real(mean(nFXYAll))],[imag(mean(nFXYAll))],...
    'ko','markerf','y')
hold on
plot([0 real(mean(nFXYAll))],...
    [0 imag(mean(nFXYAll))],'y-','linew',3)

hold off

%% Coherence spectrum
% computes a coh over (fake) windows point at
% frequency ff (previous cell) for all frequencies
% in freqvector. 

clear, clf, clc, close,
srate = 1000;
dt = 1/srate;
t = dt:dt:4;
f = 5;

clear CxySpectrum
freqvector = 0:0.1:10;
count = 0;

figure(1), clf
for ff = freqvector

count = count+1
K = exp(-1i*2*pi*ff*t);
clear nFXYAll

N = 25;

for nwindow = 1:N

phi = -deg2rad(90)+0.3*randn;
% need for noise in X and Y signals because if
% there is no noise it will appear coherence in
% other frequencies because in all frequencies the
% signal is non-existent and therefore similar and
% then there is coherence between them.
X = 3*sin(2*pi*f*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*f*t+phi)+0.3*randn(size(t));

% a single point, the point of the current window
% at the current freq 
W = hamming(length(X))';
FX = mean((W.*X).*K); 
FY = mean((W.*Y).*K);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I used this code snippet for ff = 5 and ff = 10
% to understand how the coherence is calculated
% figure(2), clf
% subplot(311)
% plot(t,X), hold on
% plot(t,Y), hold off
% title(['FX = ' num2str(mean((W.*X).*K)) '   FY = ' num2str(mean((W.*Y).*K))])
% subplot(312)
% plot(t,W.*X), hold on
% plot(t,W.*Y), hold off
% subplot(313)
% plot(t,((W.*X).*K)), hold on
% plot(t,((W.*Y).*K)), hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nFX = FX/abs(FX);
nFY = FY/abs(FY);

% coherence point of the freq ff and win nwindow
nFXY = nFX*conj(nFY);

% join the dots, one from each window = create a vector
nFXYAll(nwindow) = nFXY;

subplot(111)
    plot([real(nFXY)],[imag(nFXY)],'ko','markerf','k'); hold on
    plot([0 real(nFXY)],[0 imag(nFXY)],'k-')
    plot(exp(1i*(0:0.01:2*pi)),'k-'); hold on
    xlim([-1 1])
    ylim([-1 1])
    axis square
    title(['Nwindow = ' num2str(nwindow)])
% pause(0.005)
end

% midpoint of the coh at a freq over several windows
Cxy = abs(mean(nFXYAll));

% saves the Cxy for each freq
CxySpectrum(count) = Cxy;

    title(['Coherence = ' num2str(Cxy) '  Frequency = ' num2str(ff) ' Hz'])
    plot([real(mean(nFXYAll))],[imag(mean(nFXYAll))],...
        'ko','markerf','y'); hold on
    plot([0 real(mean(nFXYAll))],...
        [0 imag(mean(nFXYAll))],'y-','linew',3); hold off
% pause
end


%% Plot

clf
subplot(111)
    plot(freqvector,CxySpectrum); hold on
    plot(freqvector,CxySpectrum.^2); hold off
    xlabel('Frequency (Hz)')
    ylabel('Coherence')


%% Using matlab's built-in function
% mscohere % ms = magnitude squared

%%% run many times without clearing the figure and
%%% changing the Tmax to study the effect of
%%% sample size.


%setup
clear, clc, clf
srate = 1000;
dt = 1/srate;
% changing the sample size changes the coherence
% because it changes the amount of vectors I'm
% averaged over
Tmax = 32;
t = dt:dt:Tmax; 

% the added randn value makes no difference as it
% is not changing over time 
phi = -deg2rad(90)+0.5*randn;

% if we take the randn out of X and Y we will have
% coherence 1 in all frequencies this is what is
% happening in the program  
X = 3*sin(2*pi*10*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*10*t+phi)+0.3*randn(size(t));

% increasing the size of the window generates more
% precision in frequency but generates more noise
% in the others. 
windowlength = 2*srate;
overlap = 0;

% freqvector = 0:0.1:20;
% [Cxy F] = mscohere(X,Y,windowlength,overlap,freqvector,srate)
% plot(F,Cxy) 
% xlabel('Frequency (Hz)')
% ylabel('Coherence')

nfft = 2^16;
[Cxy, F] = mscohere(X,Y,windowlength,overlap,nfft,srate);
hold on
plot(F,Cxy)
xlabel('Frequency (Hz)')
ylabel('Coherence')
% legend('1', '2', '4', '8', '16', '32')
title('Coherence Spectrum')
xlim([0 20])

%% Coherogram (there is no biuld-in for this)
% "to compute a coherence I need to take my signal
% and break it into different windows, each window
% will give me a vector and the average vector
% size will give me the coherence spectrum. Only
% in this case I need to take a large window to
% compute my coherence spectrum and inside this
% window, get little windows to compute the
% vectors for the coherence spectrum" Adriano Tort

% Coherogram is like Time-Frequency Power
% Decomposition (TFD), the Cxy (PSD) is placed
% vertically and several Cxy (PSD) are calculated,
% one in each window. Here the additional
% complexity is that to compute the Cxy for a
% signal (in this case for a window) the Cxy
% requires breaking this window into smaller ones.

srate = 1000;
dt = 1/srate;
t = dt:dt:200;

phi = -deg2rad(90)+0.5*randn;
X = 3*sin(2*pi*10*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*10*t+phi)+0.3*randn(size(t));

% For each window of this 10 seconds I will
% compute the coherence spectrum of the previous
% cell. I'll compute one of these every 1 second =
% step size
win = 10*srate; % window for Cxy computation
step = 0.1*win;
Nwin = (length(t)-win)/step+1;

cohwin = 2*srate;
nfft = 2^16;
overlap = 0;

for nwin = 1:Nwin
nwin
winidx = (1:win) + (nwin-1)*step; 
[Cxy, F] = mscohere(X(winidx),Y(winidx),cohwin,overlap,nfft,srate);
% holds a coherence spectrum for each window
Coherogram(nwin,:) = Cxy;
T(nwin) = t(winidx(win/2));    
end
    
% %
clf
subplot(111)
    imagesc(T,F,Coherogram')
    axis xy
    ylim([0 20])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    h = colorbar;
    h.Label.String = 'Coherence';
    title('Coherogram')









