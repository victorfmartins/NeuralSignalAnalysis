%% Lesson 2 = 04/03/2021

%% Complex number representation in the complex plane
% creates a complex vector via Euclidean composition
% creates a complex vector via polar representation
% prints a line to the point
% prints a ball in place of the point

clf, clear, clc
x = 3;
y = 2;

% z = x + 1i*y; % replaced by z below
% polar representation
R = sqrt(x^2 + y^2);
theta = atan(y/x);
z = R*exp(1i*theta);

subplot(111)
    plot([0 real(z)],[0 imag(z)],'k-')
    hold on
    plot(real(z),imag(z),'ko','markerf','y','markers',14)
    hold off
    xlim([-5 5])
    ylim([-5 5])
    axis square
    title('Vector in Complex Plane')

%% Once Again
% sets modulus and angle directly
% print the vector again

R = 4;
theta = deg2rad(-45);
z = R*exp(1i*theta);

subplot(111)
    plot([0, real(z)],[0, imag(z)],'k-')
    hold on
    plot(real(z),imag(z),'ko','markerf','y','markers',14)
    hold off
    xlim([-5 5])
    ylim([-5 5])
    axis square
    title('Vector in Complex Plane')

%% Rotating the vector
% makes the vector rotate by varying the angle 0:2*pi
% shows vector rotating

for theta = 0:0.02:2*pi
    R = 4;
    z = R*exp(-1i*theta);

    plot([0, real(z)],[0, imag(z)],'k-')
    hold on
    plot(real(z),imag(z),'ko','markerf','y','markers',14)
    hold off
    xlim([-5 5])
    ylim([-5 5])
    axis square
    title('Vector in Complex Plane')
    pause(0.001)
end


%% Define Fourier Kernel
% sets Fourier Kernel
% shows that it is a vector rotating at a certain frequency
% shows vector rotating as time varies

% Fourier Kernel
f = 5; t = 1/1000:1/1000:1; K = exp(-1i*2*pi*f*t); % freely comment this line

f = 5;
for t=0:0.02:1
    theta = 2*pi*f*t;
    R = 4;
    z = R*exp(-1i*theta);

    plot([0, real(z)],[0, imag(z)],'k-')
    hold on
    plot(real(z),imag(z),'ko','markerf','y','markers',14)
    hold off
    xlim([-5 5])
    ylim([-5 5])
    axis square
    title(['Time = ' num2str(t) ' s'])
    pause(0.01)
end

%% Define a fake LFP (just a sinusoid)
% creates a sinusoidal LFP through the sin function and a temporal vector
% creates another LFP by joining sinusoids of different frequencies and amplitudes
% prints the LFP

% sinusoidal LFP
f = 3;
tvector = 0:0.001:1;
%LFP = sin(2*pi*f*tvector); % replaced by LFP below

% LFP composed of two sine at different frequencies
LFP = sin(2*pi*f*tvector)+0.2*sin(2*pi*(f+11)*tvector);

subplot(111)
plot(tvector,LFP)
    xlabel('Time (s)')
    ylabel('mV')
    set(gcf,'color','w')
    title('LFP')

%% Fourier transform in a single frequency
% determines the product of the signal with kernel one point at a time
% prints the signal dot and prints the product dot
% prints the average of the products - the coef. of Fourier on freq. Kernel

clf

% Fourier transform
% sum LFP(tvector)*exp(-1i*2*pi*f*tvector);

% max time
Tmax = 1;

ff = 8; % frequency of the Fourier Kernell, in Hz %was zero
f = 9; % signal frequency,in Hz

Zall = []; % this will accumulate all products

for t=0:0.02:Tmax

LFP = sin(2*pi*f*t); % had a +2

subplot(3,1,1)
    plot(t,LFP,'k.') % plotting a point. One per iteration of the for
    hold on
    xlim([0 Tmax])
    ylim([-3 3])
    title(['Time = ' num2str(t) ' s'])

subplot(3,1,[2 3])
    % Fourier Kernel
    z = exp(-1i*2*pi*ff*t); % plotting a point. One per iteration of the for
    
    % Signal x fourier kernel
    Z = LFP*z; % new dot which is the product of dot.kernel and dot.sign

    % save Z by concatenating the new product into the vector of all products
    Zall = [Zall,Z];

    % prints a line to the point and a marker in its place
    plot([0, real(Z)],[0, imag(Z)],'k-')
    hold on
    plot(real(Z),imag(Z),'ko','markerf','y','markers',7)
   
    % plotting the unitary circle
    plot(exp(1i*(0:0.001:2*pi)))

    xlim([-3 3])
    ylim([-3 3])
    axis square

    pause(0.001)
end

    % computing the fourier transform
    % FX = sum(Zall); % this is the correct definition
    FX = mean(Zall); % using the mean just for visualization

    % computing the power
    % Power = abs(FX)^2; % this and the next definition are the same
    Power = FX*conj(FX);

    % continuing on subplot(3,1,[2 3])
    plot([0, real(FX)],[0, imag(FX)],'r-')
    plot(real(FX),imag(FX),'ro','markerf','r','markers',7)
    title(['Power at frequency ' num2str(ff) ' Hz = ' num2str(Power) ] )
    hold off

subplot(3,1,1)
    hold off

%% Power spectrum = fourier at various frequencies
% runs the same code as in the previous cell, but loops over FreqVect
% saves the power of each frequency
clf
clear PowerSpectrum

f = 9;
FreqVector = 1:1:10;

for ff = FreqVector

Zall = []; % this will accumulate all products of one same Kernel

for t=0:0.01:Tmax

% sinusoidal signal at 9 Hz
LFP = sin(2*pi*f*t);
if t>0.5
    LFP = 0;
end

%%%%%%%%%%% Start of repeating the code from the previous cell
subplot(3,1,1)
    plot(t,LFP,'k.') % plotting a point. One per iteration of the for
    hold on
    xlim([0 Tmax])
    ylim([-3 3])
    title(['Time = ' num2str(t) ' s'])
    
subplot(3,1,[2 3])
    % Fourier Kernel
    z = exp(-1i*2*pi*ff*t); % is a point. One per iteration of the for
    
    % Signal x fourier kernel
    Z = LFP*z; % new dot which is the product of dot.kernel and dot.sign

    % save Z by concatenating the new product into the vector of all products
    Zall = [Zall,Z];

%     % prints a line to the point and a marker in its place
%     plot([0, real(Z)],[0, imag(Z)],'k-')
%     hold on
%     plot(real(Z),imag(Z),'ko','markerf','y','markers',7)
%    
%     % plotting the unitary circle
%     plot(exp(1i*(0:0.001:2*pi)))
% 
%     xlim([-3 3])
%     ylim([-3 3])
%     axis square
% 
%     pause(0.001)
end
 
    % computing the fourier transform
    % FX = sum(Zall); % this is the correct definition
    FX = mean(Zall); % using the mean just for visualization

    % computing the power
    % Power = abs(FX)^2; % this and the next definition are the same
    Power = FX*conj(FX);

%     % continuing subplot(3,1,[2 3])
%     plot([0, real(FX)],[0, imag(FX)],'r-')
%     plot(real(FX),imag(FX),'ro','markerf','r','markers',7)
%     title(['Power at frequency ' num2str(ff) ' Hz = ' num2str(Power) ] )
%     pause(0.5)
%     hold off
%%%%%%%%%%% End of repeating the code from the previous cell

% Save the power of each frequency
PowerSpectrum(ff)=Power;
end
   
%% Plots PowerSpectrum
% plots the PowerSpectrum produced in the previous cell

% plota o PowerSpectrum
plot(FreqVector,PowerSpectrum,'ko-','linewidth',3,...
    'markersize',14,'markerfacecolor','w')

% stem plots data conecting each point to axis
% insted of one another
stem(FreqVector,PowerSpectrum,'k-','linew',3)

xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%% end of class %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% breaking a line
clc

temp = 3+4 ...
    + 5 + 3 ...
    +2;

%% A little tour about plotting


% plotting a dot
x = 3;
y = 2;

plot(x,y,'ro','markerfacecolor','y','markersize',20)

% matlab uses RGB system for colors
plot(x,y,'ro','markerfacecolor',[0 0 0],'markersize',20)

%% Plotting more than one point

% plot with square marker
plot([1 2 3 5],[4 2 -1 2],'ysq-', ...
    'markersize',20,'markerfacecolor','m')
hold on

% plots another graph over the previous one
plot(1:5,[1:5].^2)
hold off

axis square

%% Working with subplots

clf

% subplot(Rows,Columns,PanelNumber)
subplot(2,2,1)
plot([1 2 3 5],[4 2 -1 2],'ysq-', ...
    'markersize',20,'markerfacecolor','m')

subplot(2,2,2)
plot(1:5,[1:5].^2)

subplot(2,1,2)
plot(1:5,sqrt([1:5]))

%% Declaring a vetor

vetor = [];
clear vetor2
for j=1:10
    
    vetor2(j) = j;
    vetor = [vetor,j];    
end

