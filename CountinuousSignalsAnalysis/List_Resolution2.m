%% Lista de Exercicios 2 = Aula 5 = 16/03/2021

% setup
clear, clf, clc, close, load('LFP_HG_HFO.mat')
srate = 1000; % in Hz
dt = 1/srate; Tmax = length(lfpHG)*dt; t = dt:dt:Tmax; % in s
% t = (1:length(lfpHG))*dt; % substitui os dois ultimos comandos da anterio


%% Questão 1
% 1) Compute a decomposição tempo-frequência (TFD) de cada canal, utilizando 
% janelas de 4 segundos com 50% de overlap, e com resolução numérica de 
% frequência de no mínimo 0.2 Hz.

window_length = 4*srate;
overlap = 0.5*window_length; % []
nfft = 5000; % numerical frequency resolution is srate/nfft

[~,~,~,TFDHG]=spectrogram(lfpHG,window_length,overlap,nfft,srate);
[S,F,T,TFDHFO]=spectrogram(lfpHFO,window_length,overlap,nfft,srate);


%% Questao 2
% 2) Crie dois subplots e plote em cada subplot a TFD de um canal. 
% Utilize labels adequados (‘Time (s)’, ‘Frequency (Hz)’), limite de escala Y 
% de 0 a 20 Hz, e informe como título do subplot o canal ao qual a TFD se refere.

subplot(211)
imagesc(T,F,TFDHG)

axis xy
ylim([0 20])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('lfpHG')
colorbar
caxis([0 0.02])

subplot(212)
imagesc(T,F,TFDHFO)

axis xy
ylim([0 20])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('lfpHFO')
colorbar
caxis([0 0.02])

%% Questao 3
% 3) A partir das TFDs, compute a série temporal de potência média na faixa 
% das oscilações teta (6-10 Hz). Plote essas séries temporais em diferentes 
% subplots (um pra cada canal).

Itheta = find(F > 6 & F<10);

ThetaPowerHG = mean(TFDHG(Itheta,:));
ThetaPowerHFO = mean(TFDHFO(Itheta,:));


subplot(211)
plot(T,ThetaPowerHG,'k-sq')
xlabel('Time (s)')
ylabel('\Theta Power')


subplot(212)
% hold on
plot(T,ThetaPowerHFO,'r-sq')
xlabel('Time (s)')
ylabel('\Theta Power')


%% Questao 4
% 4) Plote um scatter plot (gráfico de dispersão) utilizando os valores da 
% potência de teta de cada canal no tempo.

subplot(111)
plot(ThetaPowerHG,ThetaPowerHFO,'ko')

xlabel('\theta power electrode lfpHG')

ylabel('\theta power electrode lfpHFO')


%% Questao 5
% 5) (Opcional) Compute e plote um fitting (ajuste) linear a partir dos 
% dados do gráfico acima.

% fazer fitting linear
% y = a*x + b

x = ThetaPowerHG;
y = ThetaPowerHFO;
n = 1; % grau do polinomio, n= 1 para fitting linear

p = polyfit(x,y,n);

X = 0:0.001:0.025;
Y = polyval(p,X);

hold on
plot(X,Y,'g-','linew',3)
hold off

%% Questao 6
% 6) (Opcional) Compute a correlação de Pearson entre as séries temporais 
% das potências de teta, e informe o valor desta correlação no gráfico.
[r,pvalue] = corr(x',y');

% hold on
% text(0.25,0.0125,['Pearson'+'s r = '  num2str(r)])
% hold off
%% Questao 7
% 7) Encontre a frequência pico na banda teta de cada canal para cada janela 
% da TFD, e plote as séries temporais dessas frequências num mesmo subplot 
% (2,1,1). Em outro subplot (2,1,2), plote o histograma de distribuição das 
% frequências teta de cada canal.

% estimando a frequencia pico
plot(F,TFDHG(:,1))
[m,i] = max(TFDHG(:,1));

xlim([0 20])
ylabel('Power')
xlabel('Frequency (Hz)')
hold on
plot(F(i),m,'yo','markerf','y')
hold off

% procurando o valor maximo numa banda restrita
Itheta = find(F>=6 & F<=10);
F(Itheta)';
[m,i] = max(TFDHG(Itheta,1));
F(Itheta(i));

% computando theta peak frequency across time windows
clear ThetaFreqHG ThetaFreqHFO
for nwin = 1:size(TFDHG,2)

[m,i] = max(TFDHG(Itheta,nwin));
ThetaFreqHG(nwin) = F(Itheta(i));

[m,i] = max(TFDHFO(Itheta,nwin));
ThetaFreqHFO(nwin) = F(Itheta(i));
    
end

subplot(211)
plot(T,ThetaFreqHG,'b-*')
hold on
plot(T,ThetaFreqHFO,'r-*')
hold off

legend('lfpHG','lfpHFO')
legend('location','northeastoutside')
xlabel('Time (s)')
ylabel('Theta frequency (Hz)')


subplot(212)

histogram(ThetaFreqHG,6:0.25:10)
hold on
histogram(ThetaFreqHFO,6:0.25:10)
hold off

%%

[countHG bin] = hist(ThetaFreqHG,6:0.25:10);
[countHFO bin] = hist(ThetaFreqHFO,6:0.25:10);

bar(bin,countHG,0.3)
hold on
bar(bin+0.1,countHFO,0.3)
hold off


%% Questao 8
% 8) Repita o que foi feito nas questões 4 a 6 acima, mas, ao invés da 
% potência de teta, utilizando a série temporal das frequências de teta

subplot(111)
plot(ThetaFreqHG,ThetaFreqHFO,'ko')
xlabel('\theta frequency electrode lfpHG')
ylabel('\theta frequency electrode lfpHFO')

x = ThetaFreqHG;
y = ThetaFreqHFO;
n = 1; % grau do polinomio, n= 1 para fitting linear

p = polyfit(x,y,n);

X = 6:0.1:10;
Y = polyval(p,X);

hold on
plot(X,Y,'g-','linew',3)
hold off

[r,pvalue] = corr(ThetaFreqHG',ThetaFreqHFO')

text(8,6.5,["Pearson's r = "  num2str(r)])