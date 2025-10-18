%% Análise espectral completa (double-sided com fftshift)
clear; clc; close all; % fecha todas as figuras abertas

% --- Leitura do arquivo ---
T = readtable('EntradaFiltro.txt', 'Delimiter', ',', 'MultipleDelimsAsOne', true);
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
Y = fftshift(fft(tensao))/N;               % FFT centralizada
f = (-N/2 : N/2-1) * (fs/N);              % eixo de frequência (Hz)
magY = abs(Y);                            % magnitude normalizada
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
% Energia acumulada e f_p (90% ET) — ignorando f = 0
% ==========================================================
mid = N/2 + 1;                % índice correspondente a f = 0
f_pos = f(mid:end);           % frequências positivas
Y_pos = Y(mid:end);           % espectro positivo
E_pos = abs(Y_pos).^2;

% energia total positiva e energia DC (f=0)
E_total_pos = sum(E_pos);
E0 = E_pos(1);
perc_E0 = 100 * (E0 / E_total_pos);

% acumula energia a partir da primeira frequência > 0
E_sem0 = E_pos(2:end);
cumE = cumsum(E_sem0);
fp_idx = find(cumE >= 0.9 * sum(E_sem0), 1, 'first') + 1;
fp = f_pos(fp_idx);

fprintf('Energia em f=0: %.6e (%.2f%% da energia positiva total)\n', E0, perc_E0);

% --- Marca f_p no espectro ---
figure;
plot(f, magY, 'r', 'LineWidth', 1.2); hold on;
yl = ylim;
plot([fp fp], yl, '--k', 'LineWidth', 1.0);
grid on; xlabel('Frequência (Hz)'); ylabel('|Y(f)|');
title(sprintf('Espectro de Magnitude com f_p = %.3f Hz (90%% ET, sem f=0)', fp));
legend('Magnitude', sprintf('f_p (90%% ET) = %.3f Hz', fp));

% ==========================================================
% Exibição dos resultados
% ==========================================================
fprintf('===== RESULTADOS =====\n');
fprintf('N = %d, dt = %.12g s, fs = %.3f Hz\n', N, dt, fs);
fprintf('Energia total ET = %.6e\n', ET);
fprintf('Frequência de passagem f_p (90%% ET, sem f=0) = %.3f Hz\n\n', fp);

% ==========================================================
% Filtragem espectral e ajuste de alfa_p
% ==========================================================
janela = double(abs(f) <= fp);
Y_filtrado = Y .* janela;
x_temp = ifft(ifftshift(Y_filtrado)) * N;
alfa_p = (5/1023) / max(abs(x_temp(:)));
Y_filtrado = Y .* (alfa_p * janela);
x_filtrado = ifft(ifftshift(Y_filtrado)) * N;

fprintf('Valor calculado de alfa_p = %.6e\n', alfa_p);

figure;
plot(tempo, real(x_filtrado), 'm', 'LineWidth', 1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal no Tempo devido a um ripple constante \\alpha_p');

VR = 1 + alfa_p;
VR_dB = 20*log10(VR);
fprintf('VR = %.6f (valor linear)\n', VR);
fprintf('VR = %.6f dB\n', VR_dB);

% ==========================================================
% Ajuste de alfa_r (aplicado SÓ em |f| > fp; dentro = 0)
% ==========================================================
idx_fora = abs(f) > fp;
x_fora = ifft(ifftshift(Y .* idx_fora)) * N;
max_x_fora = max(abs(x_fora(:)));
if max_x_fora == 0
warning('Não há energia fora de ±fp. alfa_r não pode ser calculado.');
alfa_r = 0;
else
limite = 5/1023;
alfa_r = limite / max_x_fora;
end

Y_mod = zeros(size(Y));
Y_mod(idx_fora) = alfa_r * Y(idx_fora);
x_mod = ifft(ifftshift(Y_mod)) * N;

fprintf('alfa_r calculado = %.6e\n', alfa_r);
if alfa_r > 0
alfa_r_dB = 20*log10(alfa_r);
else
alfa_r_dB = -Inf;
end
fprintf('alfa_r em dB = %.3f dB\n', alfa_r_dB);
fprintf('Amplitude máxima esperada (|x_mod|_max) = %.6e\n', max(abs(x_mod(:))));

figure;
plot(tempo, real(x_mod), 'LineWidth', 1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title(sprintf('Sinal Reconstruído com \\alpha_r = %.6e (mask fora de |f|<=%.3f Hz)', alfa_r, fp));

figure;
plot(f, abs(Y_mod), 'b', 'LineWidth', 1.2); hold on;
yl = ylim;
plot([fp fp], yl, '--k', 'LineWidth', 1.0);
plot([-fp -fp], yl, '--k', 'LineWidth', 1.0);
grid on;
xlabel('Frequência (Hz)');
ylabel('|Y_{mod}(f)|');
title(sprintf('Espectro Modificado |Y_{mod}(f)| com \\alpha_r = %.6e', alfa_r));
legend('|Y_{mod}(f)|', '±f_p');
