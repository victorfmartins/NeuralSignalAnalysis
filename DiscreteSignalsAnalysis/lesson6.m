%% Lesson 6 2021/08/12 Place Cells

clc, clf, clear
load('SpkBuz.mat')

tvideo = video_dt*(1:length(X));

subplot(111)
    plot(tvideo,X); hold on
    plot(tvideo,Y); hold off
    xlabel('Time (s)')
    ylabel('Space (cm)')
    title('Animal Position')
    legend('X coordinate','Y coordinate')

%% Animation: animal position with trace

for j = 1000:10:length(X)
    % exclude unrecognised position points that were
    % set to -1
    plot(X(X>0),Y(Y>0),'b-'); hold on 
    % make the trasce: las 80 positions are marked
    plot(X(j-80:j),Y(j-80:j), 'k-', 'linew', 3)
    % show current position
    plot(X(j),Y(j),'yo','markersize',12,'markerfacecolor','y')
    hold off
    xlabel('X (cm)')
    ylabel('Y (cm)')
    title(['Time (s) = ' num2str(tvideo(j)) ' s'])
    pause(0.001)
end
set(gcf,'color','w')

%% Defining 2D space bins

% Define bin limits
binwidth = 5; % in cm
XX = 100:binwidth:260;
YY = 40:binwidth:200;

subplot(111)
    plot(X(X>0),Y(Y>0),'b-')
    set(gca,'xtick',XX,'ytick',YY)
    grid on
    axis tight
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Animal Position')

%% Computing spacial occupancy

for j=1:length(XX)-1
    
% get X values in each XX limited collum (XX is a collum boundary)
I = find(X>XX(j) & X<XX(j+1));

% Plot collum
% hold on
% plot(X(I),Y(I),'ro')
% hold off

for jj = 1:length(YY)-1

% get Y values in each YY limited collum (YY is a collum boundary)
II = find(Y>YY(jj) & Y<YY(jj+1));

% Plot row
% hold on
% plot(X(II),Y(II),'ko')
% hold off

% get retangular region of intersection of I
% colluns and II rows
III = intersect(I,II);

% get the time the animal spend in each
% retangualar bin.
occupancy(j,jj) = length(III)*video_dt;

% plot(X(X>0),Y(Y>0),'b.')
% hold on
% plot(X(III),Y(III),'yo')
% hold off
% pause(0.01)

end
end


% %

subplot(2,2,1)
    plot(X(X>0),Y(Y>0))
    xlim([100 260])
    ylim([40 200])
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Animal Position')

subplot(2,2,3)
    % +binwidth/2 so that the position of the data
    % is in the center off each square, not on its
    % upper left corner.
    % (1:end-1) to get the same dimention of occupancy
    imagesc(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,occupancy')
    caxis([0 10])
    xlim([100 260])
    ylim([40 200])
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Temporal Occupancy')
    axis xy
    % colorbar


%% Plot one red dot over animal position for every spike of cell 29

% get spike times from data of cell 29
spktimes = Raster{29}; % in seconds

% First plot rastergram for visual inspection
% of data
% subplot(2,2,2)
% plot(spktimes,ones(size(spktimes)),'k.')

% A not so elegant way of geting the position
% index at the time of each spike.
spkind = zeros(size(spktimes));
for j=1:length(spktimes) %for every spike time
% get the index of the moment in video whitch is
% closest to spike the spike time 
[mm, i]=min(abs(spktimes(j)-tvideo)); 
spkind(j) = i; % sequentially store the index
end

% A more elegant way of geting the index.
% spktimes in seconds
% video_srate in frames/second
% spktimes*video_srate give the closest frame from
% the spike activity. 
% Round to get the idex of one exact frame
% i.e. index of the point in spake where the spike happened
spkind2 = round(spktimes*video_srate);

% get unvalide points
Ix = find(X<0);
% Iy = find(Y<0) %same as above line

% get valid points from spkind by excluding from
% spikind the unvalide points
spkind  = setdiff(spkind,Ix); % not in use
spkind2 = setdiff(spkind2,Ix);

subplot(221)
hold on
plot(X(spkind2),Y(spkind2),'ro','markerf','r')
hold off

%% Of the red dots above, spacely separete them with the 2D grid for spike density plot

clear occupancy spikecount
for j=  1:length(XX)-1
% Same as before:
% get X values in each XX limited collum 
% (XX is a collum boundary)
I = find(X>XX(j) & X<XX(j+1));

% spkind2 has the index of the point in spake
% where the spike happened.
% get all the indexes of spkind2 that fall in a 2D vertical collum
Ispike = find(X(spkind2)>XX(j) & X(spkind2) < XX(j+1) );

% hold on
% plot(X(spkind2(Ispike)),Y(spkind2(Ispike)),'yo','markerf','y')
% plot(X(I),Y(I),'yo','markerf','y')
% hold off

for jj = 1:length(YY)-1
% Same as before:
II = find(Y>YY(jj) & Y<YY(jj+1));

% get all the indexes of spkind2 that fall in a 2D horizontal collum
IIspike = find(Y(spkind2)>YY(jj) & Y(spkind2) < YY(jj+1) );


% hold on
% plot(X(spkind2(IIspike)),Y(spkind2(IIspike)),'yo','markerf','y')
% plot(X(II),Y(II),'yo','markerf','y')
% hold off

% Unnecessair recomputation of occupancy
% Temporal occupancy of the animal in each retangular bin
occupancy(j,jj) = length(intersect(I,II))*video_dt;

% Number of spikes cell 29 spiked in each retangular bin
spikecount(j,jj) = length(intersect(Ispike,IIspike));

end
end

%% Replot previows plots and plot the remaining 2
clf,

h1 = subplot(221);
    plot(X(X>0),Y(Y>0)); hold on
    plot(X(spkind2),Y(spkind2),'ro','markerf','r')
    hold off
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Animal Position')

h2 = subplot(2,2,3);
    imagesc(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,occupancy')
    caxis([0 10])
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Temporal occupancy')
    axis xy

% New plot: spike counts
h3 = subplot(2,2,2);
    imagesc(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,spikecount')
    caxis([0 5])
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Spike Counts')
    axis xy

linkaxes([h1 h2 h3])
xlim([100 260])
ylim([40 200])

% % Computando o Place Field

% Since it is little relayble to calculate firing
% rate where the animal spend so few time we
% exclude those points by zeroing the spike count
% in them (or make it NaN to print white squares)
spikecount_temp = spikecount;
spikecount_temp(occupancy<1) = 0;
% % % spikecount_temp(occupancy<1) = NaN;

SpatialFiringRate = spikecount_temp./occupancy;


% % % SpatialFiringRate(isnan(SpatialFiringRate)) = -1;

% Last plot: firing rate: number of spikes over
% time in each bin.
subplot(2,2,4)
    imagesc(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,SpatialFiringRate')
    axis xy
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Firing Rate Map')
    % changing colormap to show with dots where we
    % have little data to compute firing rate
    % colorbar
    % corscale = colormap
    % corscale(1,:) = [1 1 1]
    % colormap(corscale)


%% Making firing rate plot smothier

% Contourf function alread smoth the plotby using
% a moving avarage of all the eight squares around it
subplot(2,2,4)
    contourf(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,...
        SpatialFiringRate',50,'linecolor','none')
    axis xy
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Firing Rate Map')

% a moving avarage is the same as a low pass
% filter. Instead of using it one can make a
% weighted avarage using a 2D gausian around the
% point as the weights. that is the work of imgaussfilt
% imgaussfilt is from processing images toolbox
SpatialFRSmooth = imgaussfilt(SpatialFiringRate,0.5);
SpatialFRSmooth(isnan(SpatialFRSmooth))=-1;
    imagesc(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,SpatialFRSmooth')
    colorbar
    axis xy
    xlabel('Space (cm)')
    ylabel('Space (cm)')
    title('Firing Rate Map')

% contourf(XX(1:end-1)+binwidth/2,YY(1:end-1)+binwidth/2,...
%     SpatialFRSmooth',50,'linecolor','none')


%% Showing how imgaussfilt works

M = zeros(33,33);
M(15,15) = 10;

subplot(211)
imagesc(M)
axis equal
axis tight
colorbar

subplot(212)
Msmooth = imgaussfilt(M,5); % change for 10 and for 2 to see
imagesc(Msmooth)
axis equal
axis tight
colorbar







































