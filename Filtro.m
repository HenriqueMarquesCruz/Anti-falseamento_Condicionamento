% processamento_corrente_and_filter_design.m
% Script para: FFT, energia, f_p (90%), projeto de filtro analógico K=2 (Butterworth)
% Gera figuras .png para inclusão no Overleaf.

clear; close all; clc;

%% === CONFIGURAÇÃO (edite conforme seu caso) ===
% Arquivo de dados (CSV/Excel). Pode ser CSV com uma coluna (corrente) ou duas (t, i).
dataFile = 'corrente_oscilo.csv';  
isExcel = false; % true se for .xlsx/.xls, false se .csv / .txt

% Se o arquivo NÃO contém coluna de tempo, informe fs (Hz):
fs = 100000; % taxa de amostragem hipotética (Hz) -- ajuste para seu arquivo!

% Coluna onde está a corrente (1 = primeira coluna). Se houver tempo na 1a coluna
% e corrente na 2a, ajuste colIndex = 2.
colIndex = 1;

% ADC (para sugestão de ripple/atenuação). Ajuste conforme seu DAQ
ADC_bits = 12;        % ex: 12 bits
ADC_Vref  = 3.3;      % volts

% Critérios de projeto (pode alterar)
ripple_sugerido_db = 0.5;    % alfa_p (dB) sugerido
atten_sugerido_db  = 60;     % alfa_R (dB) sugerido

% Escolha do tipo de filtro analógico (Butterworth aqui)
filter_type = 'butter'; % apenas butter implementado no script

% Parâmetros Sallen-Key: escolha C (prático) para implementação
C_choice = 10e-9; % 10 nF (você pode mudar para 4.7nF, 22nF, etc.)

% Nome base para salvar figuras
outprefix = 'fig_';

%% === LEITURA DOS DADOS ===
fprintf('Lendo dados de %s ...\n', dataFile);
if isExcel
    T = readtable(dataFile);
    data = table2array(T);
else
    data = readmatrix(dataFile);
end

% Detecta se há coluna de tempo
[nrows, ncols] = size(data);
if ncols >= 2
    t = data(:,1);
    i = data(:,colIndex);
    % Estima fs a partir do tempo (sobrescreve fs declarado)
    dt = mean(diff(t));
    fs = 1/dt;
    fprintf('Arquivo com tempo detectado. fs estimada = %.3f Hz\n', fs);
else
    i = data(:,colIndex);
    fprintf('Arquivo com apenas 1 coluna (corrente). Usando fs = %.3f Hz\n', fs);
    N = length(i);
    t = (0:N-1)'/fs;
end

% Remove média (opcional) para analisar componentes AC
i_dc = mean(i);
i_ac = i - i_dc;

%% === FFT e espectro unilaterial ===
N = length(i_ac);
T = 1/fs;
% FFT (sem normalização)
X = fft(i_ac);
% Frequências correspondentes
f = (0:N-1)*(fs/N);
% pega apenas metade positiva (one-sided)
half = 1:floor(N/2);
f_pos = f(half);
X_pos = X(half);

% Magnitude e fase
mag = abs(X_pos);
phase = angle(X_pos);

% Normalização opcional (comum é dividir por N para amplitude real)
% Se quiser amplitude em V ou A corresponde à dividir por N:
mag_norm = mag / N;

%% === Plot espectro magnitude e fase ===
figure(1); clf;
plot(f_pos, mag_norm, 'LineWidth', 1);
xlim([0 fs/2]);
xlabel('Frequência (Hz)'); ylabel('Magnitude (A) (normalizada)');
title('Espectro de magnitude da corrente');
grid on;
set(gcf,'Position',[100 100 800 400]);
saveas(gcf, [outprefix 'espectro_magnitude.png']);

figure(2); clf;
plot(f_pos, unwrap(phase), 'LineWidth', 1);
xlim([0 fs/2]);
xlabel('Frequência (Hz)'); ylabel('Fase (rad)');
title('Espectro de fase da corrente');
grid on;
set(gcf,'Position',[100 100 800 400]);
saveas(gcf, [outprefix 'espectro_fase.png']);

%% === Energia total E_T (usando X POS) ===
% Observação: energia do sinal no domínio da frequência proporcional a |X|^2.
Ecomp = abs(X_pos).^2;     % componente por bin (não normalizado)
E_T = sum(Ecomp);
fprintf('Energia total (soma |X|^2 nas frequências positivas): E_T = %.6e (unidade relativa)\n', E_T);

%% === Energia acumulada e f_p (90% de E_T) ===
cumE = cumsum(Ecomp);
target = 0.90 * E_T;
idx90 = find(cumE >= target, 1, 'first');
if isempty(idx90)
    warning('Não alcançou 90%% da energia no espectro positivo. Verifique N/normalização.');
    idx90 = length(cumE);
end
fp = f_pos(idx90);
omega_p = 2*pi*fp;
fprintf('Harmônico index h_90 = %d -> f_p = %.3f Hz (omega_p = %.3f rad/s)\n', idx90-1, fp, omega_p);

% Plota energia acumulada com linha de 90%
figure(3); clf;
plot(f_pos, cumE / E_T * 100, 'LineWidth', 1);
hold on; yline(90,'r--','90%%'); xline(fp,'k--',sprintf('f_p=%.1f Hz',fp));
xlabel('Frequência (Hz)'); ylabel('Energia acumulada (% de E_T)');
title('Energia acumulada no espectro');
grid on;
set(gcf,'Position',[100 100 800 400]);
saveas(gcf, [outprefix 'energia_acumulada.png']);

%% === Determinação de alpha_p e alpha_R baseada no ADC (sugestão) ===
LSB = ADC_Vref / (2^ADC_bits - 1);
V_noise_rms = LSB / sqrt(12);  % ruído equivalente uniform distribution approx
fprintf('ADC: %d bits, Vref=%.3f V -> LSB=%.6f V, V_noise_rms~%.6e V\n', ADC_bits, ADC_Vref, LSB, V_noise_rms);

% Sugerindo que componentes fora da banda sejam atenuadas abaixo de 3*V_noise_rms:
atten_needed_V = 3 * V_noise_rms;
% Converta para dB relativo ao pico do espectro (usa mag_norm max)
peakMag = max(mag_norm);
atten_needed_db = 20*log10(atten_needed_V / peakMag);
fprintf('Sugestão: atenuar componentes fora da banda abaixo de %.2f dB (comparado ao pico da FFT).\n', atten_needed_db);

% Use valores definidos no topo, mas mostra sugestão
fprintf('alpha_p (sugerido): %.2f dB (você configurou %.2f dB)\n', ripple_sugerido_db, ripple_sugerido_db);
fprintf('alpha_R (sugerido): %.2f dB (você configurou %.2f dB)\n', atten_sugerido_db, atten_sugerido_db);

%% === Projeto do filtro analógico (Butterworth) ordem K=2 ===
K = 2;
% Usaremos a frequência de passagem f_p obtida (pode ser substituída manualmente)
fc = fp;   % frequência de corte / passagem (Hz)
Wc = 2*pi*fc; % rad/s

% Projeto analógico Butterworth (prototype) com função butter (s-domain)
% Em MATLAB: [b,a] = butter(n, Wn, 's') retorna forma analógica
[b, a] = butter(K, Wc, 's'); % coeficientes do filtro analógico
sys = tf(b, a);
% Resposta em frequência (analógica)
w = logspace(log10(2*pi*max(1,fc/100)), log10(2*pi*max(fc*100,fs*pi)), 1000); % rad/s
[mag_sys, phase_sys] = bode(sys, w);
mag_sys = squeeze(mag_sys); % em ganho (abs)
phase_sys = squeeze(phase_sys);

% Plot resposta em dB e fase
figure(4); clf;
semilogx(w/(2*pi), 20*log10(mag_sys)); grid on;
xlabel('Frequência (Hz)'); ylabel('Ganho (dB)');
title(sprintf('Resposta do filtro analógico Butterworth K=%d (fc=%.2f Hz)', K, fc));
xlim([min(w)/(2*pi) max(w)/(2*pi)]);
saveas(gcf, [outprefix 'filtro_resposta_amplitude.png']);

figure(5); clf;
semilogx(w/(2*pi), phase_sys*180/pi); grid on;
xlabel('Frequência (Hz)'); ylabel('Fase (graus)');
title('Resposta em fase do filtro analógico');
xlim([min(w)/(2*pi) max(w)/(2*pi)]);
saveas(gcf, [outprefix 'filtro_resposta_fase.png']);

%% === Identificar f_R (frequência onde atinge alpha_R dB de atenuação) ===
% encontrar menor f tal que ganho(dB) <= -alpha_R (isto é, atenuação de alpha_R dB)
gain_db = 20*log10(mag_sys);
alpha_R = atten_sugerido_db;
idxR = find(gain_db <= -alpha_R, 1, 'first');
if ~isempty(idxR)
    fR = w(idxR)/(2*pi);
    omega_R = w(idxR);
    fprintf('Frequência de borda de rejeição f_R = %.3f Hz (atenuaçao >= %.2f dB)\n', fR, alpha_R);
else
    warning('Filtro projetado não atinge a atenuação alpha_R no range verificado.');
    fR = NaN; omega_R = NaN;
end

%% === Frequência mínima de amostragem f_S (Nyquist e fator de segurança) ===
if ~isnan(fR)
    fS_min = 2 * fR;
    fS_recom = 3 * fR; % fator de segurança 3x
    fprintf('f_S,min (Nyquist) = %.3f Hz; f_S recomendado (3x) = %.3f Hz\n', fS_min, fS_recom);
else
    fS_min = NaN; fS_recom = NaN;
end

%% === Projeto prático: Sallen-Key 2ª ordem (componentes) ===
% Para implementar o Butterworth 2ª ordem com Sallen-Key, precisamos de Q.
% Q de Butterworth 2ª ordem = 1/sqrt(2) ~ 0.7071
Q_des = 1/sqrt(2);

% Escolha de C: C1 = C2 = C_choice; R1 = R2 = R
C1 = C_choice; C2 = C_choice;
% Para Sallen-Key com ganho K, para R1=R2=R e C1=C2=C:
% Q = 1/(3 - K)  => K = 3 - 1/Q
K_sk = 3 - 1/Q_des;
if K_sk <= 1
    warning('K calculado para Sallen-Key <= 1; ajustar escolha de C/R. K_sk = %.3f', K_sk);
end
% Para desired wc = 2*pi*fc, with R1=R2=R and C1=C2=C:
% wc = 1/(R*C) => R = 1/(wc*C)
R = 1 / (Wc * C_choice);

% Implementação do ganho do amplificador (não inversor): K = 1 + Rf/Rg
% Então escolha Rg = 10k -> Rf = (K - 1) * Rg
Rg = 10e3;
Rf = (K_sk - 1) * Rg;

fprintf('\n--- Sallen-Key parameters (práticos) ---\n');
fprintf('Escolha C1 = C2 = %.3e F (%.1fnF)\n', C_choice, C_choice*1e9);
fprintf('R1 = R2 = %.3f ohm (use resistores padrão próximos)\n', R);
fprintf('Ganho SK K = %.4f => Rf = %.1f ohm, Rg = %.1f ohm (Rf/Rg = K-1)\n', K_sk, Rf, Rg);

% Exibir equação aproximada de Q obtida (checagem)
Q_check = 1/(3 - K_sk);
fprintf('Q desejado = %.4f ; Q obtido pela fórmula = %.4f\n', Q_des, Q_check);

% Gerar circuito esquemático textual (instrução)
fprintf('\nPara implementar no Sallen-Key (configuração não-inversora):\n');
fprintf('- Use op-amp em configuração não-inversora com ganho K = 1 + Rf/Rg.\n');
fprintf('- Entrada passa por R1 para nó A, nó A ligado a C1 para terra, nó A ligado a C2 para saída via R2.\n(Ver livros/ notas para topologia padrão de Sallen-Key.)\n');

%% === Salvar resultados em arquivo texto (opcional) ===
outtxt = [outprefix 'relatorio_resumo.txt'];
fid = fopen(outtxt,'w');
fprintf(fid,'Resumo do processamento e projeto\n');
fprintf(fid,'Arquivo: %s\n', dataFile);
fprintf(fid,'fs = %.3f Hz, N = %d\n', fs, N);
fprintf(fid,'E_T = %.6e (unidade relativa)\n', E_T);
fprintf(fid,'f_p (90%% energia) = %.3f Hz\n', fp);
fprintf(fid,'alpha_p (dB) usado = %.2f\n', ripple_sugerido_db);
fprintf(fid,'alpha_R (dB) usado = %.2f\n', atten_sugerido_db);
if ~isnan(fR)
    fprintf(fid,'f_R (bord. rejeição) = %.3f Hz\n', fR);
    fprintf(fid,'f_S_min (Nyquist) = %.3f Hz\n', fS_min);
    fprintf(fid,'f_S_recomendado (3x) = %.3f Hz\n', fS_recom);
else
    fprintf(fid,'f_R não determinado no range.\n');
end
fprintf(fid,'Sallen-Key: C = %.3e F, R = %.3f Ohm, K = %.3f (Rf=%.1f, Rg=%.1f)\n', C_choice, R, K_sk, Rf, Rg);
fclose(fid);
fprintf('Resumo salvo em %s\n', outtxt);

%% === Exportar figuras para Overleaf (nomes já salvos) ===
fprintf('Figuras salvas: %sespectro_magnitude.png, %sespectro_fase.png, %sfiltro_resposta_amplitude.png, %sfiltro_resposta_fase.png, %senergia_acumulada.png\n', outprefix, outprefix, outprefix, outprefix, outprefix);
