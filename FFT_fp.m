%% Análise espectral completa (double-sided com fftshift)
clear; clc; close all; % fecha todas as figuras abertas

% --- Leitura do arquivo ---
T = readtable('Osciloscopio.txt', 'Delimiter', ',', 'MultipleDelimsAsOne', true);
T = T(:, all(~ismissing(T)));
tensao = str2double(T{:,2});

% --- Vetor de tempo manual ---
tempo = -0.049280000000 : 0.000040000000 : 0.049960000000;

% Ajusta tamanhos se forem diferentes
if length(tempo) ~= length(tensao)
    warning('Tamanhos diferentes: tempo (%d) e tensao (%d). Ajustando pelo menor.', length(tempo), length(tensao));
    N = min(length(tempo), length(tensao));
    tempo = tempo(1:N);
    tensao = tensao(1:N);
end

% --- Parâmetros de amostragem ---
dt = mean(diff(tempo));
fs = 1/dt;
N = length(tensao);

% --- FFT e eixo de frequência (double-sided) ---
Y = fftshift(fft(tensao));               % FFT centralizada
f = (-N/2 : N/2-1) * (fs/N);             % eixo de frequência (Hz)
magY = abs(Y) / N;                       % magnitude normalizada
phaseY = angle(Y);

% --- Energia total ---
ET = sum(abs(Y).^2);

% --- Plot do sinal no tempo ---
figure;
plot(tempo, tensao, 'b', 'LineWidth', 1.2);
grid on; xlabel('Tempo (s)'); ylabel('Tensão (V)');
title('Sinal de Tensão vs Tempo');

% --- Espectro de magnitude e fase ---
figure;
subplot(2,1,1);
plot(f, magY, 'r', 'LineWidth', 1.2);
grid on; xlabel('Frequência (Hz)'); ylabel('|V(f)|');
title('Espectro de Magnitude (double-sided)');
subplot(2,1,2);
plot(f, unwrap(phaseY), 'k', 'LineWidth', 1.0);
grid on; xlabel('Frequência (Hz)'); ylabel('Fase (rad)');
title('Espectro de Fase (double-sided)');

% ==========================================================
% Energia acumulada e f_p (90% ET)
% ==========================================================
% energia parcial acumulando a partir de f = 0 (módulo)
mid = ceil(N/2);
E_pos = abs(Y(mid:end)).^2;           % lado positivo (0 até +fs/2)
cumE = cumsum(E_pos);
fp_idx = find(cumE >= 0.9 * sum(E_pos), 1, 'first');
fp = f(mid - 1 + fp_idx);

% --- Marca f_p no espectro ---
figure;
plot(f, magY, 'r', 'LineWidth', 1.2); hold on;
yl = ylim;
plot([fp fp], yl, '--k', 'LineWidth', 1.0);
grid on; xlabel('Frequência (Hz)'); ylabel('|Y(f)|');
title(sprintf('Espectro de Magnitude com f_p = %.3f Hz (90%% ET)', fp));
legend('Magnitude', sprintf('f_p (90%% ET) = %.3f Hz', fp));


% ==========================================================
% Exibição dos resultados
% ==========================================================
fprintf('===== RESULTADOS =====\n');
fprintf('N = %d, dt = %.12g s, fs = %.3f Hz\n', N, dt, fs);
fprintf('Energia total ET = %.6e\n', ET);
fprintf('Frequência de passagem f_p (90%% ET) = %.3f Hz\n\n', fp);