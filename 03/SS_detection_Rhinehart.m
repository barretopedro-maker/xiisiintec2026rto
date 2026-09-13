%% ========================================================================
%  IMPLEMENTAÇÃO DO MÉTODO R-STATISTIC PARA IDENTIFICAÇÃO DE ESTADO
%  ESTACIONÁRIO E TRANSITÓRIO EM PROCESSOS COM RUÍDO
%  Baseado em: Rhinehart, R.R. (2013). Automated Steady and Transient 
%  State Identification in Noisy Processes. ACC 2013.
% =========================================================================

clear; close all; clc;

%% ========================================================================
%  CASO 1: PROCESSO MONOVARIÁVEL
% =========================================================================

fprintf('\n========================================\n');
fprintf('CASO 1: PROCESSO MONOVARIÁVEL\n');
fprintf('========================================\n\n');

% Parâmetros de simulação
nSamples = 800;
dt = 1;  % intervalo de amostragem

% Geração do sinal verdadeiro com diferentes regiões
t = (1:nSamples) * dt;
true_value = zeros(1, nSamples);

% Região 1: Estado estacionário (t = 1 a 200)
true_value(1:200) = 10;

% Região 2: Rampa (t = 201 a 400)
for i = 201:400
    true_value(i) = 10 + 0.05 * (i - 200);
end

% Região 3: Estado estacionário (t = 401 a 600)
true_value(401:600) = 20;

% Região 4: Degrau (t = 601)
true_value(601:800) = 22;

% Adição de ruído Gaussiano (SNR ajustável)
rng(100); % semente
noise_std = 0.5;  % desvio padrão do ruído
measurement = true_value + noise_std * randn(1, nSamples);

% Parâmetros do filtro (recomendados pelo autor)
lambda1 = 0.2;   % fator para média móvel
lambda2 = 0.1;   % fator para variância (numerador)
lambda3 = 0.1;   % fator para variância (denominador)

% Inicialização
xf = zeros(1, nSamples);        % valor filtrado
nu2f = zeros(1, nSamples);      % variância por diferença ao filtrado
delta2f = zeros(1, nSamples);   % variância por diferença sequencial
measurement_old = 0;
R_stat = zeros(1, nSamples);    % estatística R
SS_flag = zeros(1, nSamples);   % flag de estado (1 = SS, 0 = TS)

% Valores críticos (recomendados pelo autor)
R_upper_critical = 2.5;   % limite superior (rejeita SS -> TS)
R_lower_critical = 1.0;   % limite inferior (rejeita TS -> SS)

% Loop principal
for i = 1:nSamples
    if i == 1
        % Inicialização
        xf(i) = measurement(i);
        nu2f(i) = 0;
        delta2f(i) = 0;
        measurement_old = measurement(i);
        SS_flag(i) = 1;  % assume estado estacionário inicialmente
    else
        % Equação (2): variância baseada na diferença com valor filtrado
        nu2f(i) = lambda2 * (measurement(i) - xf(i-1))^2 + (1 - lambda2) * nu2f(i-1);
        
        % Equação (1): valor filtrado
        xf(i) = lambda1 * measurement(i) + (1 - lambda1) * xf(i-1);
        
        % Equação (3): variância baseada na diferença sequencial
        delta2f(i) = lambda3 * (measurement(i) - measurement_old)^2 + (1 - lambda3) * delta2f(i-1);
        
        % Equação (4): estatística R
        R_stat(i) = (2 - lambda1) * nu2f(i) / delta2f(i);
        
        % Teste condicional (evita divisão por zero)
        if (2 - lambda1) * nu2f(i) > delta2f(i) * R_upper_critical
            SS_flag(i) = 0;  % Estado Transitório
        elseif (2 - lambda1) * nu2f(i) < delta2f(i) * R_lower_critical
            SS_flag(i) = 1;  % Estado Estacionário
        else
            SS_flag(i) = SS_flag(i-1);  % mantém estado anterior
        end
    end
    measurement_old = measurement(i);
end

% Plot dos resultados - Caso Monovariável
figure('Name', 'Caso Monovariável', 'Position', [100, 100, 1200, 800]);

subplot(3,1,1);
plot(t, true_value, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Valor Verdadeiro');
hold on;
plot(t, measurement, 'b.', 'MarkerSize', 3, 'DisplayName', 'Medição com Ruído');
plot(t, xf, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Valor Filtrado');
ylabel('Valor da Variável');
title('Processo Monovariável - Sinal e Filtro');
legend('Location', 'best');
grid on;

subplot(3,1,2);
plot(t, R_stat, 'g-', 'LineWidth', 1.5);
hold on;
yline(R_upper_critical, 'r--', 'R_{crítico superior}', 'LineWidth', 1.2);
yline(R_lower_critical, 'r--', 'R_{crítico inferior}', 'LineWidth', 1.2);
yline(1, 'k:', 'R=1', 'LineWidth', 1);
ylabel('Estatística R');
title('Evolução da Estatística R');
legend('R calculado', 'Location', 'best');
grid on;

subplot(3,1,3);
stairs(t, SS_flag, 'b-', 'LineWidth', 1.5);
ylim([-0.1, 1.1]);
ylabel('Estado');
xlabel('Tempo (amostras)');
title('Identificação de Estado (1 = SS, 0 = TS)');
grid on;
yticks([0, 1]);
yticklabels({'Transitório', 'Estacionário'});

fprintf('Análise do Caso Monovariável:\n');
fprintf('- Mínimo R: %.3f\n', min(R_stat));
fprintf('- Máximo R: %.3f\n', max(R_stat));
fprintf('- Média R em SS: %.3f\n', mean(R_stat(SS_flag==1)));

%% ========================================================================
%  CASO 2: PROCESSO MULTIVARIÁVEL (3 VARIÁVEIS)
% =========================================================================

fprintf('\n========================================\n');
fprintf('CASO 2: PROCESSO MULTIVARIÁVEL\n');
fprintf('========================================\n\n');

% Número de variáveis
nVars = 3;
nSamples_mv = 1000;

% Geração de dados para múltiplas variáveis
measurements_mv = zeros(nVars, nSamples_mv);
true_values_mv = zeros(nVars, nSamples_mv);
SS_flags_mv = zeros(nVars, nSamples_mv);
process_SS = zeros(1, nSamples_mv);

% Parâmetros para cada variável (diferentes dinâmicas)
time_constants = [0.05, 0.08, 0.03];  % taxas de variação diferentes
noise_levels = [0.3, 0.5, 0.4];        % níveis de ruído diferentes

% Geração das trajetórias
for var = 1:nVars
    % Geração do sinal verdadeiro
    true_val = zeros(1, nSamples_mv);
    
    % Regiões de transição
    true_val(1:300) = 10 + (var-1)*5;
    for i = 301:500
        true_val(i) = 10 + (var-1)*5 + time_constants(var) * (i - 300);
    end
    true_val(501:700) = 20 + (var-1)*5;
    true_val(701:1000) = 18 + (var-1)*5;
    
    true_values_mv(var, :) = true_val;
    
    % Adição de ruído (ruído independente para cada variável)
    measurements_mv(var, :) = true_val + noise_levels(var) * randn(1, nSamples_mv);
end

% Aplicação do método R-statistic para cada variável
for var = 1:nVars
    % Inicialização
    xf = zeros(1, nSamples_mv);
    nu2f = zeros(1, nSamples_mv);
    delta2f = zeros(1, nSamples_mv);
    measurement_old = measurements_mv(var, 1);
    R_stat_var = zeros(1, nSamples_mv);
    
    for i = 1:nSamples_mv
        if i == 1
            xf(i) = measurements_mv(var, i);
            nu2f(i) = 0;
            delta2f(i) = 0;
            SS_flags_mv(var, i) = 1;
        else
            nu2f(i) = lambda2 * (measurements_mv(var, i) - xf(i-1))^2 + (1 - lambda2) * nu2f(i-1);
            xf(i) = lambda1 * measurements_mv(var, i) + (1 - lambda1) * xf(i-1);
            delta2f(i) = lambda3 * (measurements_mv(var, i) - measurement_old)^2 + (1 - lambda3) * delta2f(i-1);
            R_stat_var(i) = (2 - lambda1) * nu2f(i) / delta2f(i);
            
            if (2 - lambda1) * nu2f(i) > delta2f(i) * R_upper_critical
                SS_flags_mv(var, i) = 0;
            elseif (2 - lambda1) * nu2f(i) < delta2f(i) * R_lower_critical
                SS_flags_mv(var, i) = 1;
            else
                SS_flags_mv(var, i) = SS_flags_mv(var, i-1);
            end
        end
        measurement_old = measurements_mv(var, i);
    end
    
    R_stat_mv(var, :) = R_stat_var;
end

% Equação (6): Estado do processo = produto dos estados individuais
for i = 1:nSamples_mv
    process_SS(i) = prod(SS_flags_mv(:, i));
end

% Plot dos resultados - Caso Multivariável
figure('Name', 'Caso Multivariável', 'Position', [100, 100, 1200, 1000]);

for var = 1:nVars
    subplot(nVars+1, 1, var);
    plot((1:nSamples_mv), measurements_mv(var, :), 'b.', 'MarkerSize', 2);
    hold on;
    plot((1:nSamples_mv), true_values_mv(var, :), 'k-', 'LineWidth', 1.5);
    ylabel(sprintf('Variável %d', var));
    title(sprintf('Variável %d - Medição e Valor Verdadeiro', var));
    legend('Medição com Ruído', 'Valor Verdadeiro', 'Location', 'best');
    grid on;
end

subplot(nVars+1, 1, nVars+1);
stairs((1:nSamples_mv), process_SS, 'b-', 'LineWidth', 1.5);
ylim([-0.1, 1.1]);
ylabel('Estado do Processo');
xlabel('Tempo (amostras)');
title('Estado do Processo Multivariável (1 = Todos em SS, 0 = Pelo menos um em TS)');
grid on;
yticks([0, 1]);
yticklabels({'Transitório', 'Estacionário'});

fprintf('Análise do Caso Multivariável:\n');
fprintf('- Número de variáveis: %d\n', nVars);
fprintf('- Percentual de tempo em estado estacionário: %.1f%%\n', 100 * sum(process_SS) / nSamples_mv);
fprintf('- Coordenação entre variáveis:\n');
for var = 1:nVars
    fprintf('  * Variável %d: SS em %.1f%% do tempo\n', var, 100 * sum(SS_flags_mv(var, :)) / nSamples_mv);
end

%% ========================================================================
%  FUNÇÃO AUXILIAR: MÉTODO R-STATISTIC (VERSÃO FUNÇÃO)
% =========================================================================

function [SS_flag, R_stat] = R_Statistic_Method(measurements, lambda1, lambda2, lambda3, R_upper, R_lower)
% R_Statistic_Method - Implementação do método R-statistic
%
% Entradas:
%   measurements - vetor de medições
%   lambda1, lambda2, lambda3 - fatores dos filtros (recomendado: 0.2, 0.1, 0.1)
%   R_upper, R_lower - valores críticos superior e inferior
%
% Saídas:
%   SS_flag - vetor de flags (1 = SS, 0 = TS)
%   R_stat - vetor com valores da estatística R

    n = length(measurements);
    SS_flag = zeros(1, n);
    R_stat = zeros(1, n);
    
    xf = zeros(1, n);
    nu2f = zeros(1, n);
    delta2f = zeros(1, n);
    
    measurement_old = measurements(1);
    SS_flag(1) = 1;
    
    for i = 1:n
        if i == 1
            xf(i) = measurements(i);
            nu2f(i) = 0;
            delta2f(i) = 0;
        else
            nu2f(i) = lambda2 * (measurements(i) - xf(i-1))^2 + (1 - lambda2) * nu2f(i-1);
            xf(i) = lambda1 * measurements(i) + (1 - lambda1) * xf(i-1);
            delta2f(i) = lambda3 * (measurements(i) - measurement_old)^2 + (1 - lambda3) * delta2f(i-1);
            R_stat(i) = (2 - lambda1) * nu2f(i) / delta2f(i);
            
            % Evita divisão por zero na comparação
            if delta2f(i) > 1e-10
                if (2 - lambda1) * nu2f(i) > delta2f(i) * R_upper
                    SS_flag(i) = 0;
                elseif (2 - lambda1) * nu2f(i) < delta2f(i) * R_lower
                    SS_flag(i) = 1;
                else
                    SS_flag(i) = SS_flag(i-1);
                end
            else
                SS_flag(i) = SS_flag(i-1);
            end
        end
        measurement_old = measurements(i);
    end
end

%% ========================================================================
%  ANÁLISE DE DESEMPENHO (OPCIONAL)
% =========================================================================

fprintf('\n========================================\n');
fprintf('ANÁLISE DE DESEMPENHO\n');
fprintf('========================================\n\n');

% Simulação para análise de sensibilidade
fprintf('Análise de Sensibilidade ao Ruído:\n');

noise_levels_test = [0.1, 0.5, 1.0, 2.0];
for k = 1:length(noise_levels_test)
    % Gera sinal com degrau
    true_signal = [10*ones(1, 200), 20*ones(1, 300)];
    noisy_signal = true_signal + noise_levels_test(k) * randn(1, 500);
    
    % Aplica método R-statistic
    [SS_flag_test, R_test] = R_Statistic_Method(noisy_signal, 0.2, 0.1, 0.1, 2.5, 1.0);
    
    % Tempo de detecção da transição
    transition_idx = find(SS_flag_test(201:end) == 0, 1);
    if ~isempty(transition_idx)
        detection_time = transition_idx;
    else
        detection_time = NaN;
    end
    
    fprintf('  Ruído std = %.1f: Tempo de detecção = %d amostras\n', ...
        noise_levels_test(k), detection_time);
end

fprintf('\nAnálise Concluída.\n');