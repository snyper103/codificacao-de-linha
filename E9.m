%Gabriel Alexandre de Souza Braga 
clear;
close all;
clc;
%% Prametros %%
Fs = 1e3;   % Freq. de amostragem
Ts = 1/Fs;  % Período de amostragem 

% Gera os dados binários de forma aleatória com distribuição uniforme de probabilidades
Nb = 100;   % Quantidade de bits considerados  
b = rand(1,Nb)>0.5;

% Gera o sinal unipolar NRZ já adicionando jitter
s = [];
apb = 100;  % Amostras por bit - a frequência de símbolos é Fs/apb = 1e3/100 = 10 Hz
sj = 0.01*apb;  % variância do jitter (utilizei uma distribuição normal com média 0)
for n = 1:length(b)
  if b(n) == 1
    s = [s ones(1,apb + round(sqrt(sj)*randn))];  % Adiciona simbolo '1'
  else
    s = [s zeros(1,apb + round(sqrt(sj)*randn))]; % Adiciona símbolo '0'
  end
end

t = 0:Ts:Ts*(length(s)-1);  % Tempo de sinal

%% Código de linha AMI %%
j = apb;
cl = s;
k = 1;
for i = 1:length(s)
    if s(i) == 1 && j >= apb/2
        cl(i) = k;
        j=j-1;

    elseif s(i) == 1
        cl(i) = 0;
        j=j-1;

    else
        cl(i) = 0;
        if j ~= apb
            j = apb;
            k = -k;
        end
    end
    
    if j == -1
        j = apb;
        k = -k;
    end
end

% Passa o sinal codificado em linha pelo canal h(t)
wc = 20;                          % Freq. de corte de 20 Hz
[bf,af] = butter (10, wc/(Fs/2)); % Canal como um FPB butterworth
x0 = filter(bf,af,cl);              

% Adiciona ruído AWGN
snr = 30; % Relação sinal/ruído em [dB]
x = awgn(x0, snr);

Tapb = -0.5:1/199:0.5;  % tempo para o diagrama de olho

%% Plota gráficos no domínio do tempo %%
% Cria figura
figure1 = figure('PaperOrientation', 'landscape', 'PaperUnits', 'centimeters',...
    'PaperType', 'A4',...
    'WindowState', 'maximized',...
    'Color', [1 1 1],...
    'Renderer', 'painters');

% Cria subplot dos dados binários
subplot1 = subplot(3, 1, 1, 'Parent', figure1);
hold(subplot1, 'on');

% Cria plot
stem(b, 'DisplayName', 'Dados binários', 'Parent', subplot1, 'LineWidth', 3,...
    'Color', [0.00,0.45,0.74]);

% Cria rotulo y e x
ylabel('b[n]', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('n', 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% Cria titulo
title('Dados binários');

% Define limites do plot, para x e y
xlim(subplot1, [0 Nb]);
ylim(subplot1, [-0.1 1.1]);

% Liga as grades e etc
box(subplot1, 'on');
grid(subplot1, 'on');
hold(subplot1, 'off');

% Define as propriedades restantes dos eixos
set(subplot1, 'AlphaScale', 'log', 'ColorScale', 'log', 'FontName',...
    'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold', 'GridAlpha', 0.5,...
    'LineWidth', 1.5, 'MinorGridAlpha', 0.5);

% Plota da codificado unipolar NRZ e AMI
% Cria subplot
subplot2 = subplot(3, 1, 2, 'Parent', figure1);
hold(subplot2, 'on');

% Cria plot
plot(t, s, 'DisplayName', 'NRZ', 'Parent', subplot2, 'LineWidth', 3,...
    'Color', [0.00,0.45,0.74]);
plot(t, cl, 'DisplayName', 'AMI', 'Parent', subplot2, 'LineWidth', 3,...
    'Color', [0.85 0.33 0.098], 'LineStyle', '-.');

% Cria rotulo y e x
ylabel('s(t)', 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% Cria titulo
title('Codificado unipolar NRZ e bipolar AMI');

% Define limites do plot, para x e y
xlim(subplot2, [0 Ts*(length(s)-1)]);
ylim(subplot2, [-1.25 1.25]);

% Liga as grades e etc
box(subplot2, 'on');
grid(subplot2, 'on');
hold(subplot2, 'off');

% Define as propriedades restantes dos eixos
set(subplot2, 'AlphaScale', 'log', 'ColorScale', 'log', 'FontName',...
    'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold', 'GridAlpha', 0.5,...
    'LineWidth', 1.5, 'MinorGridAlpha', 0.5);

% Define legendas
legend1 = legend(subplot2,'show');
set(legend1,'Location','southeast','LineWidth',1,'FontSize',16);

% Plota do sinal pelo canal (FPB + AWGN)
% Cria subplot
subplot3 = subplot(3, 1, 3, 'Parent', figure1);
hold(subplot3, 'on');

% Cria plot
plot(t, x, 'DisplayName', '(FPB + AWGN)', 'Parent', subplot3, 'LineWidth', 3,...
    'Color', [0.00,0.45,0.74]);

% Cria rotulo y e x
ylabel('x(t)', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('t', 'FontWeight', 'bold',...
    'FontName', 'Times New Roman');

% Cria titulo
title('Após passar pelo canal (FPB + AWGN)');

% Define limites do plot, para x e y
xlim(subplot3, [0 Ts*(length(s)-1)]);
ylim(subplot3, [-1.5 1.5]);

% Liga as grades e etc
box(subplot3, 'on');
grid(subplot3, 'on');
hold(subplot3, 'off');

% Define as propriedades restantes dos eixos
set(subplot3, 'AlphaScale', 'log', 'ColorScale', 'log', 'FontName',...
    'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold', 'GridAlpha', 0.5,...
    'LineWidth', 1.5, 'MinorGridAlpha', 0.5);

%% Diagrama de olho %%
eyediagram(x,2*apb,1,apb); % Utilize esta função pronta somente para testar
% Cria figura
figure2 = figure('PaperOrientation', 'landscape', 'PaperUnits', 'centimeters',...
    'PaperType', 'A4',...
    'WindowState', 'maximized',...
    'Color', [1 1 1],...
    'Renderer', 'painters');

% Cria subplot da magnitude espectral do sinal modulado
axes1 = axes('Parent', figure2);
hold(axes1, 'on');

plot2 = zeros(1, round(length(x)/(2*apb),0));
% Cria plot
for i = 1:length(x)/(2*apb)
    plot2(i) = plot(Tapb, x((i-1)*2*apb+1:i*2*apb), 'LineWidth', 3, 'Color', [0.00,0.45,0.74]);
end

% Cria rotulo y e x
ylabel('Amplitude', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Tempo', 'FontWeight', 'bold',...
    'FontName', 'Times New Roman');

% Cria titulo
title('Diagrama de olho');

% Define limites do plot, para x e y
xlim(axes1, [-0.5 0.5]);
ylim(axes1, [-1.5 1.5]);

% Liga as grades e etc
box(axes1, 'on');
grid(axes1, 'on');
hold(axes1, 'off');

% Define as propriedades restantes dos eixos
set(axes1, 'AlphaScale', 'log', 'ColorScale', 'log', 'FontName',...
    'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold', 'GridAlpha', 0.5,...
    'LineWidth', 1.5, 'MinorGridAlpha', 0.5);
