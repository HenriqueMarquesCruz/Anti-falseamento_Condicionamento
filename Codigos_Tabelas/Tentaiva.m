%% Análise espectral completa (double-sided com fftshift)
clear; clc; close all;

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
else
    N = length(tensao);
end

% --- Força tudo a ser vetor coluna ---
tensao = tensao(:);
tempo = tempo(:);

% --- Parâmetros de amostragem ---
dt = mean(diff(tempo));
fs = 1/dt;
fprintf('Amostras: %d | Fs = %.2f Hz\n', N, fs);

% --- FFT e eixo de frequência alinhados com fftshift ---
Y = fft(tensao);
Ys = fftshift(Y);

if mod(N,2) == 0
    f = (-N/2 : N/2 - 1) * (fs / N);
else
    f = (-(N-1)/2 : (N-1)/2) * (fs / N);
end
f = f(:); % garante vetor coluna

% --- Cálculo da energia espectral ---
E_total = sum(abs(Ys).^2);
[~, idx0] = min(abs(f));  % índice de f=0
E0 = abs(Ys(idx0)).^2;
E_pos = sum(abs(Ys(idx0:end)).^2);

fprintf("Energia total: %.3e | Energia em f=0: %.3e (%.2f%% da energia positiva)\n", ...
        E_total, E0, 100*E0/E_pos);

% --- Determinação de f_p (90%% da energia positiva, sem DC) ---
E_meta = 0.9 * (E_pos - E0);
acum = cumsum(abs(Ys(idx0+1:end)).^2);  % ignora DC
idx_rel = find(acum >= E_meta, 1);
idx_fp = idx0 + idx_rel;
fp = f(idx_fp);

fprintf('Frequência de passagem f_p = %.3f Hz\n', fp);

% --- Máscaras de passagem e rejeição ---
janela = double(abs(f) <= fp);
idx_fora = double(abs(f) > fp);

fprintf('Soma(janela)+Soma(idx_fora) = %d (esperado = N=%d)\n', sum(janela)+sum(idx_fora), N);

% --- Filtragem espectral ---
Y_inband = Ys .* janela;
Y_outband = Ys .* idx_fora;

% --- Transformadas inversas ---
x_in = real(ifft(ifftshift(Y_inband)));
x_out = real(ifft(ifftshift(Y_outband)));

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

fprintf('\nVR = %.6f (valor linear)\n', VR);
fprintf('VR = %.6f dB\n', VR_dB);
fprintf('alfa_r calculado = %.9e\n', alfa_r);
fprintf('alfa_r em dB = %.3f dB\n', alfa_r_dB);

fprintf('\nmax_x_in = %.6e, max_x_out = %.6e\n', max(abs(x_in)), max(abs(x_out)));
fprintf('rms_in = %.6e, rms_out = %.6e\n', rms(x_in), rms(x_out));
fprintf('alfa_p = %.12e, alfa_r = %.12e\n', alfa_p, alfa_r);

if abs(alfa_p - alfa_r)/max(alfa_p, alfa_r) < 0.05
    warning('alfa_p e alfa_r estão numericamente próximos — verifique máscaras e sinais reconstruídos.');
else
    fprintf('OK: alfa_p e alfa_r distintos.\n');
end

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
