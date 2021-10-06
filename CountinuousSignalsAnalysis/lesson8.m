%% Lesson 8 = 30/03/2021; Convolution

%% 'same' convolution and 'default' convolution
% graphic display of signal, kernel and conv 
clear

X = [0 2 1 1 -1 7 -2 -3 8 1 -1 2 -3 4];
K = [1 1 -2];

% Using build in function
ConvXK = conv(X,K);
ConvXKB = conv(X,K,'same'); %yellow part has the same length as signal

subplot(311)
plot(2:15,X,'bo-')
xlim([0 17])

subplot(312)
plot(K,'ro-')
xlim([0 17])

subplot(313)
plot(ConvXK,'ko-')
hold on
plot(2:15,ConvXKB,'yo-')
hold off
xlim([0 17])

%% Application example: moving average
% moving average is taking the sum of a sequence
% of n values and dividing by n and doing this by
% moving the points always adding one at the
% leading end of the sequence and taking one at
% the trailing end of the sequence

clf

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

% random signal (noise)
LFP = randn(size(t));

subplot(311)
plot(t,LFP)
xlabel('time (s)')

% ordem = #points in kernel
ordem = 25;
% moving average
% Constant Kernel has zero centered half gaussian 
% frequency spectrum. Therefore a low filter (low
% pass)   
K = ones(1,ordem)/ordem;

% Convol = conv(LFP,K);
Convol = conv(LFP,K,'same');

hold on
plot(t,Convol,'k-','linew',3)
hold off
legend('signal', 'conv')

subplot(312)
[Pxx, F] = pwelch(LFP,length(LFP),[],2^16,srate);
hold off
plot(F,Pxx)
xlim([0 500])
ylim([0 0.02])
xlabel('freq (hz)')


subplot(313)
[Pcc, F] = pwelch(Convol,length(Convol),[],2^16,srate);
hold off
plot(F,Pcc)
xlim([0 500])
ylim([0 0.02])
xlabel('freq (hz)')



%% Power Spectrum
% Power Spectrum with the same window length

% the larger the size of the K, the more localized is its power 
K = ones(1,50);

subplot(311)
[Pxx, F] = pwelch(LFP,length(K),[],2^16,srate);
plot(F,Pxx)
xlim([0 500])
xlabel('freq (hz)')
ylabel('Pxx')

subplot(312)
[Pkk, F] = pwelch(K,length(K),[],2^16,srate);
plot(F,Pkk)
xlim([0 500])
xlabel('freq (hz)')
ylabel('Pkk')

subplot(313)
Pcc = Pxx.*Pkk;
plot(F,Pcc)
xlim([0 500])
xlabel('freq (hz)')
ylabel('Pcc')


%% Shows the edge effect and its size depending on the order 

% clf
srate = 1000; dt = 1/srate; t = dt:dt:2;

% LFP = randn(size(t));
LFP(500:1500) = LFP(500:1500) + 30*exp(-t(500:1500));

% THE SIZE OF THE EDGE EFFECT DEPENDS ON THE SIZE OF THE KERNELL
% does kernel size influence angle shift? Why 90 degrees? 
ordem = 200;
K = ones(1,ordem)/ordem;

Convol = conv(LFP,K,'same');

subplot(111)
plot(t,LFP)
hold on
plot(t,Convol,'r-','linew',3)
plot([0.2 0.2],[-5 20],'r--','linew',2)
hold off

xlabel('time (s)')

%% Powers: Pxx, Pkk, Pxx.*Pkk
% conv in time equals multiplication in frequency
%%% Why does the filter work? P q ff = 10 filter f = 40? 
clf, clear, clc

% setup
srate = 1000;
dt = 1/srate;
t = dt:dt:1;

LFP = sin(2*pi*10*t)+1.5*sin(2*pi*5*t)+0.5*sin(2*pi*40*t);
LFP = LFP + randn(size(t));

% Changing the order from 200 to 500 shows the multiplication effect 
order = 100;
K = sin(2*pi*40*t(1:order));
norm = sum(K.^2); % to avoid adding or removing energy

subplot(311)
[Pxx, F] = pwelch(LFP,length(LFP),[],2^16,srate);
plot(F,Pxx)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Pxx')

subplot(312)
[Pkk, F] = pwelch(K/norm,length(K),[],2^16,srate);
plot(F,Pkk)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Pkk')

subplot(313)
plot(F,Pkk.*Pxx)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Power convol')


%% Show Kernel walking on signal generating convolution result

figure(1)
LFP= sin(2*pi*10*t);
LFP = LFP + 0.3*randn(size(t));
% padding with zeros
LFP = [zeros(1,order) LFP zeros(1,order)];
t = (1:length(LFP))*dt;

for j = 0:2:length(LFP)-order
    
subplot(211)
plot(t,LFP)
hold on
plot(t((1:order)+j),K(order:-1:1),'r-','linew',2)
hold off
ylim([-4 4])

subplot(212)

plot(t(order/2+j),sum(K(order:-1:1).*LFP((1:order)+j))/norm,'ko')
hold on
xlim([0 t(end)])
ylim([- 4 4])
pause(0.001)
end

%% Eliminating convolution phase shifts
% convolution causes 90o phase shift conv why 90o
% phase shift conv? 
% why do you have to normalize by sum(k.^2) by
% energy? how much does the order influence?  

srate = 1000; dt = 1/srate; t = dt:dt:2;
LFP = sin(2*pi*8*t)+0.5*sin(2*pi*20*t);

order = 500;
K = sin(2*pi*8*t(1:order));
norm = sum(K.^2); % to avoid adding or removing energy

Convol = conv(LFP,K/norm,'same');

% to avoid phase shifts, convolve the 
% convolution with the inverted kernell
Convol2 = conv(Convol,K(end:-1:1)/norm,'same');

subplot(411)
hold on
plot(t,LFP,'r-')
hold off

[Pxx, ~] = pwelch(LFP, length(K), [], 2^14, srate);
[Pkk, F] = pwelch(K, length(K), [], 2^14, srate);
subplot(412)
plot(F,Pxx/sum(Pxx))
hold on
plot(F,Pkk/sum(Pkk))
hold off
xlim([0 25])

subplot(4,1,[3 4])
plot(t, LFP,'k-')
hold on
plot(t,Convol,'r-')
plot(t,Convol2,'g-')
hold off

%% Studying the 90th shift
% The shift happens because when the tip of the K
% starts to leave the zero pad section the
% convolution point (in the middle of the K)
% starts to raise order/2 before (pi)    

clf, clear
srate = 1000; dt = 1/srate; t = dt:dt:0.5;

order = 100;
K = sin(2*pi*10*t(1:order));
norm = sum(K.^2); % to avoid adding or removing energy

LFP= sin(2*pi*10*t);;
% padding with zeros
LFP = [zeros(1,order) LFP zeros(1,order)];
t = (1:length(LFP))*dt;


figure(1)
for j = 0:2:length(LFP)-order
    
subplot(211)
plot(t,LFP)
hold on
plot(t((1:order)+j),K(order:-1:1),'k-','linew',2)
plot(t(order/2+j),sum(K(order:-1:1).*LFP((1:order)+j))/norm, 'ko')
% write code that indicates the point-to-point
% value of the product K(.)*LFP(.)  
hold off
ylim([-4 4])

subplot(212)

plot(t(order/2+j),sum(K(order:-1:1).*LFP((1:order)+j))/norm, 'ko')
hold on
xlim([0 t(end)])
ylim([- 4 4])
pause
end

%% stain the area between 2 functions matlab

clear, clf
x = linspace(-2,2);
curve1 = x.^3;
curve2 = x;
plot(x, curve1, 'r', 'LineWidth', 2);
hold on;
plot(x, curve2, 'b', 'LineWidth', 2);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g');






















