% =========================================================================
% CÓDIGO UNIFICADO: FILTRAGEM + DETECÇÃO DE ESTADO ESTACIONÁRIO (RHINEHART)
% =========================================================================
clear; close all; clc;

% 1. PARÂMETROS
n_turbines = 100; 
SP_ref_MW = 0.31;
nsim = 500; 

% Fator do filtro de ruído EWMA principal
alpha = 0.15; 

% Parâmetros do Diagnóstico de Rhinehart (Baseado em xf_alpha = 0.2)
alpha1 = 0.2; 
alpha2 = 0.2; 
alpha3 = 0.2;
R_threshold = 3.0; % Limiar típico para transição (Abaixo disso = Steady State)

% Inicialização das variáveis
P_bruto = zeros(1, nsim); 
P_limpo = zeros(1, nsim); 
SS_flag = ones(1, nsim);
R_stat = zeros(1, nsim); 

% Variáveis internas de Rhinehart
xf = zeros(1, nsim);
v1 = zeros(1, nsim);
v2 = zeros(1, nsim);

% 2. LOOP PRINCIPAL (FILTRAGEM + DIAGNÓSTICO)
for k = 1:nsim
    % A. Simulação (Sinal com Ruído Intencional)
    % Adicionando um degrau artificial no meio (k=250) para testar a detecção de Rhinehart
    if k < 250
        P_bruto(k) = 31000 + 1500 * randn(); 
    else
        P_bruto(k) = 35000 + 1500 * randn(); % Mudança de estado (Transient State)
    end
    
    % B. FILTRO DE RUÍDO PRINCIPAL (EWMA)
    if k == 1
        P_limpo(k) = P_bruto(k);
        xf(k) = P_limpo(k);
    else
        P_limpo(k) = alpha * P_bruto(k) + (1 - alpha) * P_limpo(k-1);
    end
    
    % C. IMPLEMENTAÇÃO COMPLETA DO DIAGNÓSTICO DE RHINEHART
    if k > 1
        % 1. Média móvel filtrada (Cálculo de xf que já estava no seu código)
        xf(k) = alpha1 * P_limpo(k) + (1 - alpha1) * xf(k-1);
        
        % 2. Variância baseada na diferença da média (v1)
        v1(k) = alpha2 * (P_limpo(k) - xf(k))^2 + (1 - alpha2) * v1(k-1);
        
        % 3. Variância baseada em diferenças consecutivas (v2)
        v2(k) = alpha3 * (P_limpo(k) - P_limpo(k-1))^2 + (1 - alpha3) * v2(k-1);
        
        % 4. Cálculo da Estatística R (R-Statistic)
        if v2(k) ~= 0
            R_stat(k) = (2 - alpha1) * v1(k) / v2(k);
        else
            R_stat(k) = 1;
        end
        
        % 5. Diagnóstico Final: 1 = Steady State (SS), 0 = Transient State (TS)
        if R_stat(k) > R_threshold
            SS_flag(k) = 0; % Estado Transiente
        else
            SS_flag(k) = 1; % Estado Estacionário
        end
    else
        R_stat(k) = 1;
        SS_flag(k) = 1;
    end
end

% =========================================================================
% 3. VISUALIZAÇÃO DAS FIGURAS SOLICITADAS
% =========================================================================

% --- FIGURA 1 ---
figure('Color','w');
plot(P_bruto, 'r', 'LineWidth', 0.5);
title('Potência Ativa (Sinal c/ Ruído)');
xlabel('Tempo'); ylabel('Potência (kW)');
grid on;

% --- FIGURA 2 ---
figure('Color','w');
plot(xf, 'b', 'LineWidth', 2);
title('Média Filtrada (x\_f)');
xlabel('Tempo'); ylabel('Valor de x_f');
grid on;

% --- FIGURA 3 ---
figure('Color','w');
plot(R_stat, 'k', 'LineWidth', 1.5);
title('R\_stat, ''k''');
xlabel('Tempo'); ylabel('R-Value');
grid on;

% --- FIGURA 4 ---
figure('Color','w');
plot(R_stat, 'm', 'LineWidth', 1.5);
hold on;
yline(R_threshold, 'r--', 'Limiar Critico', 'LineWidth', 1.5);
title('Estatística R-Rhinehart');
xlabel('Tempo'); ylabel('Valor R');
legend('R-Statistic', 'Limiar');
grid on;

% --- FIGURA 5 (Parte 1: Diagnóstico) ---
figure('Color','w');
plot(SS_flag, 'g', 'LineWidth', 2);
ylim([-0.2 1.2]);
title('Diagnóstico Final: 1=SS, 0=TS');
xlabel('Tempo'); ylabel('Estado');
yticks([0 1]); yticklabels({'0 (Transiente - TS)', '1 (Estacionário - SS)'});
grid on;

% --- FIGURA 6 (Conforme pedido no final da sua lista) ---
figure('Color','w');
plot(P_bruto, 'k', 'LineWidth', 0.5, 'DisplayName', 'Potência Bruta (c/ Ruído)');
hold on;
plot(P_limpo, 'b', 'LineWidth', 2, 'DisplayName', 'Potência Limpa (Filtrada)');
yline(31000, 'r--', 'Alvo Inicial', 'LineWidth', 1.5);
title('Processamento de Sinal: Remoção de Ruído para Despacho');
xlabel('Tempo'); ylabel('Potência (kW)');
legend('Location', 'best'); 
grid on;