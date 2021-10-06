%% Lesson 2 = 2021/07/27; Spike Clustering
% Recomended material: http://www.scholarpedia.org/article/Spike_sorting
% Open problens that this class will comment about: 
% spike detection, feature extraction and clustering
clear, clc, clf
load('rawdata3.mat')
srate = 20000; % in Hz
dt = 1/srate;
t = dt*(1:length(data)); % in seconds

% Arbitrary treshold (see last class)
threshold = mean(data)+2*std(data);

% Class 1 was concerned about building a findpeaks
% function.
[~,Ispikes2]= findpeaks(data,'MinPeakHeight',threshold);

% pré alocation of memory for claimed speed
% improvment
waveshape = zeros(length(Ispikes2),41);
tic
for j=1:length(Ispikes2)
    waveshape(j,:) = data(Ispikes2(j)-20:Ispikes2(j)+20);
end
toc
subplot(211)
    plot(t,data','k')
    axis tight
    title('MUA')
    ylabel('mV')
    xlabel('Time (s)')
    
subplot(212)
    plot(-20:20,waveshape','k')
    axis tight
    title('Waveshape')
    ylabel('mV')
    xlabel('Lags from spike peak')
    
%% Feature choice
% For clustering one need unicly determaining data
% to separete spikes unequivocaly. 
% One can choose wave value at a time when each
% neuron has a distinct value from others. 

% one could chose to work with all wave poins in a
% multidimentional space I guess, but the
% reducionality of dimentions here is to chose
% only 2 points to simplify understanding. In
% lecture 3 a better method is introduced:
% dimentionality reduction by PCA

subplot(2,2,[1,2])
plot(-20:20,waveshape','k')
    axis tight
    title('Waveshape')
    ylabel('mV')
    xlabel('Lags from spike peak')

% Feature0 is not a good one for clustering because
% it daoes not distinguish clearly between neurons
% as one can see in its histogram
Feature0 = waveshape(:,5);
Feature1 = waveshape(:,17);
Feature2 = waveshape(:,35);
Feature3 = waveshape(:,25);

subplot(223)
    % hist(Feature0,-2:0.02:1)
    hist(Feature1,-2:0.02:1)
    
subplot(224)
    plot(Feature1,Feature2,'ko')
    % plot3(Feature1,Feature2,Feature3,'ko')
    xlabel('Feature 1')
    ylabel('Feature 2')
    % zlabel('Feature 3')

%% Clustering via k-means method
% 1)Define a number of clusters
% 2)Give a random vector with dim=#clusters to each
% cluster. They define the center of the cluster.
% 3)Compute the norm from all points to each center.
% 4)The center closest to a point is the cluster for
% that point.
% 5)Atualize the center of the cluster by the mean
% point of all the point of that cluster.
% 6)Repete previows 3 steps until the diference of
% the maximum center displacement is below some
% delta.

% Arbitrarely define centers of dim=3 for 3
% clusters (random values)
Center1 = [0,0];
Center2 = [0.01,0];
Center3 = [-0.01,0];

% Do only 5 repetitions of steps 3 to 5
for j=1:5

% get the distence of each point to each center
temp = [repmat(Center1,length(Feature1),1) - [Feature1,Feature2]];
dist1 = sqrt(temp(:,1).^2+temp(:,2).^2);

temp = [repmat(Center2,length(Feature1),1) - [Feature1,Feature2]];
dist2 = sqrt(temp(:,1).^2+temp(:,2).^2);

temp = [repmat(Center3,length(Feature1),1) - [Feature1,Feature2]];
dist3 = sqrt(temp(:,1).^2+temp(:,2).^2);

% get to witch center the distence is the smallest
[m, idx] = min([dist1,dist2,dist3]');

% replot all points labeling them by color. Each
% point get the collor of it's closest center.
plot(Feature1,Feature2,'ko')
hold on
plot(Feature1(idx==1),Feature2(idx==1),'go')
plot(Feature1(idx==2),Feature2(idx==2),'bo')
plot(Feature1(idx==3),Feature2(idx==3),'ro')

% Mark current position of the centers
plot(Center1(1),Center1(2),'kx','markersize',20)
plot(Center2(1),Center2(2),'kx','markersize',20)
plot(Center3(1),Center3(2),'kx','markersize',20)
hold off

% Make the new center the center os the points of
% its cluster
Center1n = [mean(Feature1(idx==1)),mean(Feature2(idx==1))];
Center2n = [mean(Feature1(idx==2)),mean(Feature2(idx==2))];
Center3n = [mean(Feature1(idx==3)),mean(Feature2(idx==3))];

% Get the change in center position from
% iteraction to iteraction
DeltaCenter1(j) = norm(Center1n-Center1);
DeltaCenter2(j) = norm(Center2n-Center2);
DeltaCenter3(j) = norm(Center3n-Center3);

Center1 = Center1n;
Center2 = Center2n;
Center3 = Center3n;

if max([DeltaCenter1(j),DeltaCenter2(j),DeltaCenter3(j)])<0.0001
   disp('convergiu')
    break
end
subplot(224)
    title(['interaction ' num2str(j)])
    xlabel('Feature 1')
    ylabel('Feature 2')
% pause

end
subplot(224)
    title('Clusters')
    xlabel('Feature 1')
    ylabel('Feature 2')
%% Usin built in for K-means
% One of the field problems is to decide how many
% centers to declare at start.

clf
% if 3<numk<6 uncoment apropriate lines bellow
numk = 3;
[idx, C] = kmeans([Feature1,Feature2],numk);

subplot(111)
    plot(Feature1,Feature2,'ko')
    hold on
    plot(Feature1(idx==1),Feature2(idx==1),'go')
    plot(Feature1(idx==2),Feature2(idx==2),'bo')
    plot(Feature1(idx==3),Feature2(idx==3),'ro')
    plot(C(1,1),C(1,2),'kx','markersize',20)
    plot(C(2,1),C(2,2),'kx','markersize',20)
    plot(C(3,1),C(3,2),'kx','markersize',20)

    % plot(Feature1(idx==4),Feature2(idx==4),'yo')
    % plot(Feature1(idx==5),Feature2(idx==5),'mo')
    % plot(Feature1(idx==6),Feature2(idx==6),'co')
    % plot(C(4,1),C(4,2),'kx','markersize',20)
    % plot(C(5,1),C(5,2),'kx','markersize',20)
    % plot(C(6,1),C(6,2),'kx','markersize',20)

    hold off
    title('Clusters')
    xlabel('Feature 1')
    ylabel('Feature 2')

%% Plot waveshape of each cluster and its mean

clf
subplot(111)
    plot(-20:20,waveshape(idx==1,:)','g-')
    hold on
    plot(-20:20,waveshape(idx==2,:)','b-')
    plot(-20:20,waveshape(idx==3,:)','r-')

    plot(-20:20,mean(waveshape(idx==1,:))','k-','linew',3)
    plot(-20:20,mean(waveshape(idx==2,:))','k-','linew',3)
    plot(-20:20,mean(waveshape(idx==3,:))','k-','linew',3)

    % plot(-20:20,waveshape(idx==4,:)','y-')
    % plot(-20:20,waveshape(idx==5,:)','m-')
    % plot(-20:20,waveshape(idx==6,:)','c-')
    % 
    % plot(-20:20,mean(waveshape(idx==4,:))','k-','linew',3)
    % plot(-20:20,mean(waveshape(idx==5,:))','k-','linew',3)
    % plot(-20:20,mean(waveshape(idx==6,:))','k-','linew',3)

    hold off
    axis tight
    title('Waveshape')
    ylabel('mV')
    xlabel('Lags from spike peak')

%% Manual Clustering

clf
% draw by clicking around a visual cluster in the
% plot a closed poligon. Each point inside the
% poligon will belong to the same cluster

% first clustering
plot(Feature1,Feature2,'ko')
title('Escolha Cluster 1')
[x,y] = getline(gca,'close');
Idx1 = inpolygon(Feature1,Feature2,x,y);
title('Cluster 1')
hold on
plot(x,y,'g-')
plot(Feature1(Idx1),Feature2(Idx1),'go')
hold off

% Second clustering
title('Escolha Cluster 2')
[x,y] = getline(gca,'close');
Idx2 = inpolygon(Feature1,Feature2,x,y);
title('Cluster 2')
hold on
plot(x,y,'b-')
plot(Feature1(Idx2),Feature2(Idx2),'bo')
hold off

% Third clustering
title('Escolha Cluster 3')
[x,y] = getline(gca,'close');
Idx3 = inpolygon(Feature1,Feature2,x,y);
title('Cluster 3')
hold on
plot(x,y,'r-')
plot(Feature1(Idx3),Feature2(Idx3),'ro')
hold off

%% and ploting waveshape of clustered spikes

plot(waveshape(Idx1,:)','g-'), hold on
plot(waveshape(Idx2,:)','b-')
plot(waveshape(Idx3,:)','r-')

plot(mean(waveshape(Idx1,:))','k-','linew',3)
plot(mean(waveshape(Idx2,:))','k-','linew',3)
plot(mean(waveshape(Idx3,:))','k-','linew',3), hold off

%% Plotando o rastergrama
% time position of spikes of each cluster


h1 = subplot(211);
plot(t(Ispikes2(Idx1)),ones(1,sum(Idx1)),'gv'), hold on
plot(t(Ispikes2(Idx2)),ones(1,sum(Idx2))+1,'bv')
plot(t(Ispikes2(Idx3)),ones(1,sum(Idx3))+2,'rv'), hold off
ylim([0 4])
ylabel('Neuron #')
xlabel('Time (s)')
set(gca,'ytick',1:3)

h2 = subplot(212);
plot(t,data), hold on
plot(t(Ispikes2(Idx1)),data(Ispikes2(Idx1)),'gv')
plot(t(Ispikes2(Idx2)),data(Ispikes2(Idx2)),'bv')
plot(t(Ispikes2(Idx3)),data(Ispikes2(Idx3)),'rv'), hold off

linkaxes([h1 h2],'x')







