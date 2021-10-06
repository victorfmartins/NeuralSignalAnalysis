%% Lesson 3 = 2021/07/29; PCA
% o autovetor da matriz de covariancia vai aontar
% para a direção de maior dvariancia possivel
% Vari?ncia

clear, clc, clf
% Revising Varience and Covarience
X = rand(1,6);
Y = [2,4,0,0,4,2];

mean(X);
mean(Y);

var(X);
var(Y);
% var(Y,1);

% Explicit formula for variance and covariance
variance = sum((Y-mean(Y)).^2)/(length(Y)-1);
covariance = sum((Y-mean(Y)).*(X-mean(X)))/(length(Y)-1);

% Covariance matrix
% elements representation: aij. 
% a11 is var(X)
% a22 is var(Y)
% a12 and a21 is the covariance of X and Y
C = cov(X,Y); % generate a matrix


%% Similarity between Covariance and Corralation

X = randn(1,100);
Y = 4*X + randn(1,100);

[r, p] = corr(X',Y');
cov(X,Y)/(std(X)*std(Y));


%% Resume class 2b

clear, clc
X = randn(1,100);
Y = 4*X + randn(1,100);
Vetores = [X;Y];

Direcao = [1 4];
Direcao = Direcao/norm(Direcao);


plot(X,Y,'ko')
hold on
plot(3*[-Direcao(1) Direcao(1)],...
    3*[-Direcao(2) Direcao(2)])
hold off

ProjAll = Direcao*Vetores;
ProjvectAll = Direcao'*ProjAll;
% ProjvectAll = ProjAll'*Direcao;

hold on
plot(ProjvectAll(1,:),ProjvectAll(2,:),'go')
hold off

title(['Variance = ' num2str(var(ProjAll))])

% axis equal
xlim([-5 5])


%% Use PCA as a redution of dimentionality technique for waveshapes

clear, clc, clf
X = randn(1,100);
Y = 4*X + randn(1,100);
Vetores = [X;Y];

% computa matriz de covariancia
C = cov(X,Y);
% PC = Principal Components
% PC are the eigenvectors of the cov matrix
% lambda are the eigenvalues
% this eigenvalues are the variance. The first aii
% is lowest variance and the last aii the highest.
% We call the highest var the PC1 and the lowest
% the last PC
[PC,lambda] = eig(C);

% Sort from lower to higher values to match PC indexes
[variancias, I]=sort(diag(lambda),'descend');

PC1 = PC(:,I(1)); % can use same index after sorting
PC2 = PC(:,I(2));

subplot(3,2,[1 4])
    plot(X,Y,'ko'); hold on
    plot(3*[-PC1(1) PC1(1)],3*[-PC1(2) PC1(2)],'linew',2)
    plot(3*[-PC2(1) PC2(1)],3*[-PC2(2) PC2(2)],'linew',2); hold off
    axis equal

Proj1 = PC1'*Vetores;
Proj2 = PC2'*Vetores;

var(Proj1) % = variancias(1)
var(Proj2) % = variancias(2)

subplot(3,2,5)
    plot(Proj1,Proj2,'ko')
    xlabel('PC1')
    ylabel('PC2')
    axis equal

% scree shows the % of variance in each PC
scree = variancias/sum(variancias);

subplot(3,2,6)
    plot(100*scree,'b-o','markerfacecolor','b')
    xlabel('Component Number')
    ylabel('% Variance')
    set(gca,'xtick',1:length(PC))
    xlim([0 length(PC)+1])

    
%% applicando PCA para classificar neuronios

clear, clc, clf

% Toy model. N for neuron. codes waveshape of 4 points
N1 = randn(100,4); % N spikes vs waveshape 
N2 = randn(100,4); % 100 spikes, only 4 pints each 

% to not let all random make them diverge
N1(:,2) = N1(:,2)+5;
N2(:,2) = N2(:,2)-5;
N1(:,3) = N1(:,3)-5;
N2(:,3) = N2(:,3)+5;

M = [N1;N2];

C = cov(M);
[PC,lambda] = eig(C);
% Sort from lower to higher values to match PC indexes
[variancias, I]=sort(diag(lambda),'descend');
PC1 = PC(:,I(1)); % can use same index after sorting
PC2 = PC(:,I(2));
% PC3 = PC(:,I(3));
% PC4 = PC(:,I(4));


subplot(3,2,[1 2])
    plot(M','k')
subplot(323)
    stem(PC1)
    % xlim([0 3])
    title('PC1')
subplot(324)
    stem(PC2)
    % xlim([0 3])
    title('PC2')

    
% % Seeing the clusters
Feature1 = M*PC1;
Feature2 = M*PC2;


subplot(3,2,5)
    plot(Feature1,Feature2,'ko')
    xlabel('PC1')
    ylabel('PC2')


% % Scree
scree = variancias/sum(variancias);


subplot(3,2,6)
    plot(100*scree,'b-o','markerfacecolor','b')
    xlabel('Component Number')
    ylabel('% Variance')
    set(gca,'xtick',1:length(PC))
    xlim([0 length(PC)+1])


%% Same as before, but with more waveshape points and neurons

N1 = randn(100,20); % N spikes vs waveshape 
N2 = randn(100,20); 
N3 = randn(100,20);
N4 = randn(100,20);

N1(:,8:12) = N1(:,8:12)+5;
N2(:,8:12) = N2(:,8:12)-5;
N3(:,8:12) = N3(:,8:12)-5;
N3(:,4:6)  = N3(:,4:6)+5;
N4(:,10:12) = N4(:,10:12)+5;

M = [N1;N2;N3;N4];

C = cov(M);
[PC,lambda] = eig(C);
[variancias, I]=sort(diag(lambda),'descend');
PC1 = PC(:,I(1));
PC2 = PC(:,I(2));
PC3 = PC(:,I(3));
PC4 = PC(:,I(4));


subplot(3,2,[1 2])
    plot(M','k')
    axis tight
    title('Fake Waveshape')
    ylabel('mV')
subplot(323)
    stem(PC1)
    % xlim([0 3])
    title('PC1')
subplot(324)
    stem(PC2)
    % xlim([0 3])
    title('PC2')

    
% % Seeing the clusters
Feature1 = M*PC1;
Feature2 = M*PC2;
Feature3 = M*PC3;


subplot(325)
    plot(Feature1,Feature2,'ko')
    xlabel('PC1')
    ylabel('PC2')
    title('Clusters')
%     plot(Feature1,Feature3,'ko') % or
%     xlabel('PC1')
%     ylabel('PC3')
%     plot3(Feature1,Feature2,Feature3,'ko') % or
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')


% % Scree
scree = variancias/sum(variancias);


subplot(326)
    plot(100*scree,'b-o','markerfacecolor','b')
    xlabel('Component Number')
    ylabel('% Variance')
    set(gca,'xtick',1:length(PC))
    xlim([0 length(PC)+1])
    title('Scree')

    
%% Code from aula2
% reintroducing the real waveshape

clear, clc, clf
load('rawdata3.mat')
srate = 20000; % in Hz
dt = 1/srate;
t = dt*(1:length(data));

% Arbitrary treshold (see last class)
threshold = mean(data)+2*std(data);

% Class 1 was concerned about building a findpeaks
% function.
[temp,Ispikes2]= findpeaks(data,'MinPeakHeight',threshold);

% pré alocation of memory for claimed speed
% improvment
waveshape = zeros(length(Ispikes2),41);
tic
for j=1:length(Ispikes2)
    waveshape(j,:) = data(Ispikes2(j)-20:Ispikes2(j)+20);
end
toc
subplot(111)
    plot(-20:20,waveshape','k')
    axis tight
    title('Waveshape')
    ylabel('mV')
    xlabel('Lags from spike peak')
%%% Here is the starting point of this class %%%


%% Repeating this class code now on the real waveshape

[PC,lambda] = eig(cov(waveshape));

% Sort from lower to higher values to match PC indexes
[variancias, I]=sort(diag(lambda),'descend');
PC1 = PC(:,I(1));
PC2 = PC(:,I(2));

% % Seen PC1 and PC2
subplot(3,2,[1 2])
    plot(-20:20,waveshape','k')
    axis tight
    title('Waveshape')
    ylabel('mV')
    xlabel('Lags from spike peak')
subplot(323)
    stem(PC1)
    % xlim([0 3])
    title('PC1')
subplot(324)
    stem(PC2)
    % xlim([0 3])
    title('PC2')


% % Seeing the clusters
Feature1 = waveshape*PC1;
Feature2 = waveshape*PC2;


subplot(325)
    plot(Feature1,Feature2,'ko')
    xlabel('PC1')
    ylabel('PC2')
    title('Clusters')
    

% % Scree
scree = variancias/sum(variancias);


subplot(326)
    plot(100*scree,'b-o','markerfacecolor','b')
    xlabel('Component Number')
    ylabel('% Variance')
    set(gca,'xtick',1:length(PC))
    xlim([0 length(PC)+1])
    title('Scree')
% realize that only PC1 and PC2 acount for almost
% 90% of the varience in the waveshapes


% From here there is the k-means clustering
% algoritm seen in aula2. Thus PCA works as a
% dimentionality reduction technique.







