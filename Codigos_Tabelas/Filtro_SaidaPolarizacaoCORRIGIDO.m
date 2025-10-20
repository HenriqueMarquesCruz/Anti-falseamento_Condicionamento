%% Análise espectral completa (double-sided com fftshift)
clear; clc; close all;

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
else
    N = length(tensao);
end

% --- Força tudo a ser vetor coluna ---
tensao = tensao(:);
tempo = tempo(:);

% --- Parâmetros de amostragem ---
dt = mean(diff(tempo));
fs = 1/dt;
fprintf('Amostras: %d | fs = %.2f Hz\n', N, fs);

% --- FFT e eixo de frequência alinhados com fftshift ---
Y = fft(tensao)/N;
Ys = fftshift(Y);

if mod(N,2) == 0
    f = (-N/2 : N/2 - 1) * (fs / N);
else
    f = (-(N-1)/2 : (N-1)/2) * (fs / N);
end
f = f(:); % garante vetor coluna

% --- Cálculo da energia espectral ---
E_total = sum(abs(Ys).^2);
[~, idx0] = min(abs(f)); 
E0 = abs(Ys(idx0)).^2;
E_pos = sum(abs(Ys(idx0:end)).^2);

fprintf("Energia total: %.3e | Energia em f=0: %.3e (%.2f%% da energia positiva)\n", ...
        E_total, E0, 100*E0/E_pos);

% --- Plot do sinal no tempo ---
figure;
plot(tempo, tensao, 'b', 'LineWidth', 1.2);
grid on; xlabel('Tempo (s)'); ylabel('Tensão (V)');
title('Sinal de Tensão vs Tempo');

% --- Magnitude e fase ---
magY = abs(Ys);
phaseY = angle(Ys);

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

% --- Determinação de f_p (90%% da energia positiva, sem DC) ---
E_meta = 0.9 * (E_pos - E0);
acum = cumsum(abs(Ys(idx0+1:end)).^2);  % ignora DC
idx_rel = find(acum >= E_meta, 1);
idx_fp = idx0 + idx_rel;
fp = f(idx_fp);

% --- Marca f_p no espectro ---
figure;
plot(f, magY, 'r', 'LineWidth', 1.2); hold on;
yl = ylim;
plot([fp fp], yl, '--k', 'LineWidth', 1.0);
grid on; xlabel('Frequência (Hz)'); ylabel('|Y(f)|');
title(sprintf('Espectro de Magnitude com f_p = %.3f Hz (90%% ET)', fp));
legend('Magnitude', sprintf('f_p (90%% ET) = %.3f Hz', fp));

fprintf('Frequência de passagem f_p = %.3f Hz\n', fp);

% --- Máscaras de passagem e rejeição ---
janela = double(abs(f) <= fp);
idx_fora = double(abs(f) > fp);


% --- Filtragem espectral ---
Y_inband = Ys .* janela;
Y_outband = Ys .* idx_fora;

% --- Transformadas inversas ---
x_in = real(ifft(ifftshift(Y_inband) * N));
x_out = real(ifft(ifftshift(Y_outband) * N));

% --- Força a forma de vetor coluna (para evitar erro) ---
x_in = x_in(:);
x_out = x_out(:);

% --- Cálculo dos fatores alfa_p e alfa_r ---
alfa_p = (1/1023) / max(abs(x_in(:)));
alfa_r = (1/1023) / max(abs(x_out(:)));

% --- Resultados ---
VR = 1 + alfa_r;
VR_dB = 20*log10(VR);
alfa_r_dB = 20*log10(alfa_r);

fprintf('VR = %.6f dB\n', VR_dB);
fprintf('alfa_r calculado = %.9e\n', alfa_r);


fprintf('alfa_p = %.12e\n', alfa_p);



% --- Gráficos de verificação ---
figure;
subplot(3,1,1);
plot(f, abs(Ys)); grid on;
xlabel('Frequência (Hz)'); ylabel('|Y(f)|');
title('Espectro original');

subplot(3,1,2);
plot(f, janela, 'b', f, idx_fora, 'r');
xlabel('Frequência (Hz)'); ylabel('Máscaras');
legend('Passa-baixa', 'Rejeição', 'Location', 'best');
title('Máscaras espectrais');

subplot(3,1,3);
plot(f, abs(Y_inband), 'b'); hold on;
plot(f, abs(Y_outband), 'r');
xlabel('Frequência (Hz)'); ylabel('|Y_{mod}(f)|');
legend('Banda de passagem', 'Banda rejeitada');
title('Espectros modificados');
grid on;