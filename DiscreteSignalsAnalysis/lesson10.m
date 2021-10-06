%% Lesson 10 = 2021/09/09 ROC Curve
%% Concept: points in the ROC Curve

% set random data once
clear, clc, clf
N_trials = 100;
FR_Carro = 2*randn(1,N_trials) + 8; % firing rate for car  stimulus
FR_Face = 2*randn(1,N_trials) + 12; % firing rate for face stimulus

% get the current state of the random generators.
% every random number generate after a rng(s)
% calling will output the same random values
s = rng; 
for limiar = 0:20;

TP = sum(FR_Face>=limiar); % absolute number of true  positives
FN = sum(FR_Face< limiar); % absolute number of false negatives

TPR = TP/(FN+TP); % true positive rate (aka sensitivity)

FP = sum(FR_Carro>=limiar); % false positive
TN = sum(FR_Carro< limiar); % true negative

FPR = FP/(FP+TN); % false positive rate (aka 1-specificity)

subplot(121)
    % restart the random generator to state s (of
    % before the loop so that the plots dont
    % change in each iteration  
    rng(s);
    % all data lies in a vertical line and many
    % markers can overlap due to the proximity of
    % points. To make them more visible assign a
    % horizontal displacement to each one using
    % randn function
    plot(2*ones(size(FR_Carro))+ ...
      0.1*randn(size(FR_Carro)),FR_Carro,'rv'); hold on
    plot(1*ones(size(FR_Face))+ ...
      0.1*randn(size(FR_Face)), FR_Face, 'kv')
    plot([0 3],[limiar limiar ],'m--','linew',3); hold off
    ylim([0 20])
    ylabel('Firing Rate (Hz)')
    xlabel('Stimulus')
    set(gca,'xtick',1:2,'xticklabel',{'Face';'Car'},'fontsize',12)
    xlim([0 3])

subplot(122)
    plot([0 1],[0 1],'k--'); hold on
    plot(FPR,TPR,'ko','markerface','k','markersize',12);
    xlabel('FPR (1-Specificity)')
    ylabel('TPR (Sensitivity)')
    axis square
    set(gca,'fontsize',12)

pause(0.1)
end
    hold off


%% ROC Curve

% get all points, put them in order and use each
% one as treshhold for fine grain precision
AllFR = [FR_Carro,FR_Face];
threshold = sort(AllFR);

clear TPR FPR
counts = 0;
for limiar = threshold
    counts = counts+1;
    
TP = sum(FR_Face>=limiar); % absolute number of true  positives
FN = sum(FR_Face< limiar); % absolute number of false negatives

TPR(counts) = TP/(FN+TP); % true positive rate (aka sensitivity)

FP = sum(FR_Carro>=limiar); % false positive
TN = sum(FR_Carro< limiar); % true negative

FPR(counts) = FP/(FP+TN); % false positive rate (aka 1-specificity)   
end

% AUC = Area Under the Curve
% function trapz makes a trapesoidal (or midd point) integration
AUC = trapz(FPR(end:-1:1),TPR(end:-1:1));

subplot(122)
    hold on
    plot(FPR,TPR,'linew',3)
    hold off
    title(['ROC-AUC = ' num2str(AUC)])









