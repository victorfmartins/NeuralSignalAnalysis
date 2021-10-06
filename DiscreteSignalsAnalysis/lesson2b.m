%% Lesson 2b = 2021/07/27; PCA - Principle Component Analysis

% Linear algebra revision
% For some matrixes a given vector will not be
% rotated (those are eigenvectors) only compressed
% of expended. The compretion/expantion amount is
% the eigenvalue
clear

% Matrix
A = [2 1
    1 2];

% Vectors
v = [1 ; 1];
w = A*v;

plot([0 v(1)],[0 v(2)],'k-','linew',2)
hold on
plot([0 w(1)],[0 w(2)],'r-','linew',2)
hold off

xlim([-3 3])
ylim(xlim())

axis square
% axis equal

%% Polar representation of vectors

A = [1 2
    2 3];

% Each point represents a vector
% The circle helps to see what hapens to vectors
% in all angles

Circulo = exp(1i*(0:0.01:2*pi));
C = [real(Circulo);imag(Circulo)];

xlim([-3 3])
ylim(xlim())

% C Transfor by A
T = A*C;

% Below are eigenvectors for matriz A
% ponto=102; % it is a vector that only expands (do not rotate)
ponto=260; % it is a vector that only shrinks (do not rotate)
    
% Print the regular circle
plot(C(1,:),C(2,:));
hold on

% print one line as a vector instead of point
plot([0 C(1,ponto)],[0 C(2,ponto)],'k-')
plot([0 T(1,ponto)],[0 T(2,ponto)],'r-')

% Prints the transformed circle
plot(T(1,:),T(2,:),'r')
hold off

% xlim([-3 3])
% ylim(xlim())
% axis square
axis equal

% Theorem: The area of the transformed object is
% equal to que area of the original object times
% the determinant of the transformation matriz

% In our case here, det(A) = -1 so the area does
% not changes:
polyarea(C(1,:),C(2,:));
polyarea(T(1,:),T(2,:));

% Computing eigenvalues and eingenvectors
[V,lambda]= eig(A);
v1 = V(:,1);
v2 = V(:,2);
lambda = diag(lambda);

%% Projection of V over W (vectors)

clear, clf, clc
% setup two vectors
V = [1 3];
W = [1 1];

W = W/norm(W); % gets the direction of W
Proj = V*W'; % gets the projection of V over W
% tranpose W only for matrix multiplication reasons

plot([0 V(1)],[0 V(2)],'k-','linew',2)
hold on
% plot([0 W(1)],[0 W(2)],'r-','linew',2)

plot(5*[-W(1) W(1)],5*[-W(2) W(2)],'r-','linew',2)
plot(Proj*[0 W(1)],Proj*[0 W(2)],'b-','linew',2)
hold off

xlim([-3 3])
ylim(xlim())

% axis square
axis equal

%% Projection of many vectors

clear, clf, clc
% creates 100 random and noise [X 3X] points
X = randn(1,100);
Y = 3*X + randn(1,100);
plot(X,Y,'ko')
hold on

% Chosen direction to witch every point will be
% projected
Direcao = [1 1];
Direcao = Direcao/norm(Direcao);
plot(5*[-Direcao(1) Direcao(1)],5*[-Direcao(2) Direcao(2)],'r-','linew',2)

% Vectorize the points X and Y and take one vector
Vetores = [X;Y];
v = Vetores(:,1);

% get magnitude of projectionand vectorizes it
Proj = Direcao*v;
Projvect = Proj*Direcao;

% do the same to all 100 vectors
ProjAll = Direcao*Vetores;
ProjvectAll = Direcao'*ProjAll;

% plot projections in green
plot(ProjvectAll(1,:),ProjvectAll(2,:),'go')
plot(v(1),v(2),'ko','markerfacecolor','k')
plot(Projvect(1),Projvect(2),'ko','markerfacecolor','b')
hold off
axis equal














