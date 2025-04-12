close all;
clear;
clc;

load('data.mat');

Ts = 0.1; % Temps d'échantillonnage défini dans l'énoncé

figure;
subplot(2,1,1); % Premier graphique pour l'entrée u
plot(time, u);
title('Signal d''entrée u(t)');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2); % Second graphique pour la sortie y
plot(time, y);
title('Signal de sortie y(t)');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;



% Analyse de la première montée (t=0 à t=100s)
idx_start = 1; % L'échelon commence au premier point (t=0)
time_start = time(idx_start);

% Trouver l'index où u redescend pour la première fois (fin du premier plateau)
% Cela correspond à t=100s
idx_end_plateau = find(diff(u) < -0.5, 1); 
if isempty(idx_end_plateau)
    idx_end_plateau = length(u); % Sécurité si jamais il n'y a qu'une montée
end

% Valeurs de l'entrée u
u_step_val = u(idx_start); % Valeur de l'échelon u (devrait être 1)
u_prev_val = 0; % On suppose que u était 0 avant t=0
delta_u = u_step_val - u_prev_val; % Variation de l'entrée (devrait être 1)

% Valeurs de la sortie y
y_init = y(idx_start); % Valeur initiale de y (proche de 0)

% Estimer y_final en prenant la moyenne sur la fin du plateau pour lisser le bruit
% Par exemple, sur les points entre t=80s et t=100s (indices approx 801 à 1001)
idx_avg_start = max(idx_start + 1, idx_end_plateau - round(20/Ts)); % Début moyennage (ex: 20s avant la fin)
idx_avg_end = idx_end_plateau;                                     % Fin moyennage
y_final_est = mean(y(idx_avg_start : idx_avg_end)); 
delta_y = y_final_est - y_init; % Variation de la sortie

% 1. Estimation du Gain Statique K
K_est_graph = delta_y / delta_u;

% 2. Estimation de la Constante de Temps tau
% Calculer la valeur cible à 63%
y_63_target = y_init + 0.63 * delta_y;

% Trouver l'index où y dépasse y_63_target pour la première fois
% On cherche seulement dans la zone de montée (idx_start à idx_end_plateau)
idx_63 = find(y(idx_start:idx_end_plateau) >= y_63_target, 1) + idx_start - 1;

% Calculer tau si l'index est trouvé
if isempty(idx_63)
    warning('N''a pas pu trouver le point y à 63%% sur la première montée.');
    tau_est_graph = NaN; % Mettre NaN si non trouvé
    time_63 = NaN;
else
    time_63 = time(idx_63);
    tau_est_graph = time_63 - time_start; % Calcul de tau
end

% 3. Calcul des paramètres discrets a1, b1 (modèle ARX(1,1,1))
if ~isnan(tau_est_graph) && tau_est_graph > 1e-9 % Vérifier que tau est valide et positif
    a1_est_graph = -exp(-Ts/tau_est_graph);
    b1_est_graph = K_est_graph * (1 + a1_est_graph);
else
    a1_est_graph = NaN;
    b1_est_graph = NaN;
    warning('Tau non valide, a1 et b1 non calculés.');
end

% Affichage des résultats de l'estimation graphique
fprintf('--- Estimation Graphique (1er ordre - basée sur Figure 1) ---\n');
fprintf('Analyse de la première montée (t=%.1fs à t=%.1fs)\n', time_start, time(idx_end_plateau));
fprintf('  Delta u = %.3f\n', delta_u);
fprintf('  y_init = %.3f, y_final_est = %.3f, Delta y = %.3f\n', y_init, y_final_est, delta_y);
fprintf('  Gain Statique Estimé K = %.3f\n', K_est_graph);
fprintf('  Valeur cible y @ 63%% = %.3f\n', y_63_target);
if ~isnan(tau_est_graph)
   fprintf('  Temps atteint @ 63%% = %.3f s (Index %d)\n', time_63, idx_63);
   fprintf('  Constante de temps Estimée tau = %.3f s\n', tau_est_graph);
   fprintf('  Paramètre Estimé a1 = %.4f\n', a1_est_graph);
   fprintf('  Paramètre Estimé b1 = %.4f\n', b1_est_graph);
else
   fprintf('  Constante de temps tau: Non trouvée ou invalide\n');
   fprintf('  Paramètres a1, b1: Non calculés\n');
end
fprintf('------------------------------------------------------------\n');

N = length(time);
N_est = floor(N / 2); % Utiliser la première moitié pour l'estimation
N_val = N - N_est;   % Utiliser la seconde moitié pour la validation

% Données d'estimation
u_est = u(1:N_est);
y_est = y(1:N_est);
time_est = time(1:N_est);

% Données de validation
u_val = u(N_est+1:N);
y_val = y(N_est+1:N);
time_val = time(N_est+1:N);

% --- Identification ARX(1,1,1) par Moindres Carrés ---
na = 1; % Ordre de A(z)
nb = 1; % Ordre de B(z)
nk = 1; % Retard d (u(k-nk))

% Construire la matrice Phi (lignes = phi(k)') et le vecteur Y pour l'estimation
Phi_matrix = [];
Y_vector = [];
start_idx = max(na, nb+nk-1)+1; % Commence à k=2 pour ARX(1,1,1)
for k = start_idx : N_est
    phi_k_T = []; % Régresseur ligne phi(k)'
    % Termes AR (-y(k-1)...-y(k-na))
    for i = 1:na
        phi_k_T = [phi_k_T, -y_est(k-i)];
    end
    % Termes X (u(k-nk)...u(k-nb-nk+1))
    for i = 0:nb-1
         phi_k_T = [phi_k_T, u_est(k-nk-i)];
    end
    Phi_matrix = [Phi_matrix; phi_k_T]; % Ajoute la ligne phi(k)'
    Y_vector = [Y_vector; y_est(k)];
end

% Calculer theta_hat_1 (LS) : theta = (Phi' * Phi)^-1 * Phi' * Y
try
    theta_hat_1 = (Phi_matrix' * Phi_matrix) \ (Phi_matrix' * Y_vector); 
    
    fprintf('\n--- Résultats Étape 8 : Identification ARX(1,1,1) par LS ---\n');
    fprintf('Paramètres estimés par LS:\n');
    fprintf('  a1 = %.6f\n', theta_hat_1(1));
    fprintf('  b1 = %.6f\n', theta_hat_1(2));
    fprintf('Paramètres estimés graphiquement (rappel):\n');
    fprintf('  a1 = %.6f\n', a1_est_graph); % Utilise la valeur de l'étape 4
    fprintf('  b1 = %.6f\n', b1_est_graph); % Utilise la valeur de l'étape 4

    % --- Simulation sur données de validation ---
    y_sim_1 = zeros(size(y_val));
    
    % Utilisation de l'approche 'predict' (utilise y passés réels pour prédire y(k))
    % C'est plus simple à implémenter ici et souvent utilisé pour le Fit
    Y_val_vector = [];
    Phi_val_matrix = [];
    for k = start_idx : N_val % Attention aux indices relatifs à u_val/y_val
         % Indices globaux dans u/y pour récupérer les valeurs passées
         idx_glob = N_est + k; 
         phi_k_T_val = [];
         % Termes AR
         for i = 1:na
             phi_k_T_val = [phi_k_T_val, -y(idx_glob-i)]; % y réel passé
         end
         % Termes X
         for i = 0:nb-1
              phi_k_T_val = [phi_k_T_val, u(idx_glob-nk-i)]; % u réel passé
         end
         Phi_val_matrix = [Phi_val_matrix; phi_k_T_val];
         Y_val_vector = [Y_val_vector; y_val(k)]; % y_val(k) = y(idx_glob)
    end
    
    % Prédictions sur l'ensemble de validation
    y_sim_1 = Phi_val_matrix * theta_hat_1; 

    % Comparaison graphique
    figure;
    plot(time_val(start_idx:end), Y_val_vector, 'b', time_val(start_idx:end), y_sim_1, 'r--');
    title('Validation du Modèle ARX(1,1,1)');
    xlabel('Temps (s)');
    ylabel('Sortie y');
    legend('Sortie Mesurée (Validation)', 'Sortie Prédite (Modèle 1er Ordre)');
    grid on;
    xlim([time_val(start_idx) time_val(end)]); % Ajuster les limites pour la partie validée

    % Calcul du Fit (%)
    Error_1 = Y_val_vector - y_sim_1;
    Fit_1 = 100 * (1 - norm(Error_1) / norm(Y_val_vector - mean(Y_val_vector)));
    fprintf('Fit du modèle ARX(1,1,1) sur données de validation: %.2f %%\n', Fit_1);
    fprintf('------------------------------------------------------------\n');

catch ME
    fprintf('\nErreur lors du calcul LS pour ARX(1,1,1): %s\n', ME.message);
    fprintf('Cela peut arriver si la matrice Phi''*Phi est mal conditionnée.\n');
    % Arrêter ici si l'identification échoue
    return; 
end

% --- Etape 9 : Calcul de la fonction de coût ---

% Rappel: theta_hat_1, Phi_matrix, Y_vector ont été calculés à l'étape 8
% sur les données d'estimation.

% Calcul du coût pour theta_hat_1 (le minimum trouvé par LS)
Y_pred_est_1 = Phi_matrix * theta_hat_1; % Prédictions sur l'ensemble d'estimation
Error_est_1 = Y_vector - Y_pred_est_1;
Cost_LS_1 = sum(Error_est_1.^2);

fprintf('\n--- Résultats Étape 9 : Fonction de Coût ARX(1,1,1) ---\n');
fprintf('Coût J(theta_hat_1) sur données estimation: %.6f\n', Cost_LS_1);

% Calcul du coût pour 3 thetas aléatoires
rng('default'); % Pour reproductibilité
fprintf('Coûts pour des thetas aléatoires (proches de theta_hat_1):\n');
for i = 1:3
    % Générer theta aléatoire proche de theta_hat_1 (ex: +/- 25% par composante)
    theta_rand = theta_hat_1 .* (1 + 0.5*(rand(size(theta_hat_1)) - 0.5)); 
    
    Y_pred_rand = Phi_matrix * theta_rand;
    Error_rand = Y_vector - Y_pred_rand;
    Cost_rand = sum(Error_rand.^2);
    fprintf('  Coût J(theta_rand_%d = [%.4f; %.4f]): %.6f\n', i, theta_rand(1), theta_rand(2), Cost_rand);
end

% Optionnel: Calculer le coût pour les paramètres graphiques
theta_graph = [a1_est_graph; b1_est_graph]; % Issu de l'étape 4
if ~any(isnan(theta_graph))
    Y_pred_graph = Phi_matrix * theta_graph;
    Error_graph = Y_vector - Y_pred_graph;
    Cost_graph = sum(Error_graph.^2);
    fprintf('Coût J(theta_graphique = [%.4f; %.4f]): %.6f\n', theta_graph(1), theta_graph(2), Cost_graph);
end
fprintf('------------------------------------------------------------\n');

% --- Etape 10 : Identification ARX(2,2,1) par Moindres Carrés ---
na = 2; % Ordre de A(z)
nb = 2; % Ordre de B(z)
nk = 1; % Retard d

% Reconstruire Phi_matrix et Y_vector pour ARX(2,2,1)
% Utiliser les mêmes données d'estimation (u_est, y_est)
Phi_matrix_2 = [];
Y_vector_2 = [];
start_idx = max(na, nb+nk-1)+1; % Commence à k=3 pour ARX(2,2,1)
for k = start_idx : N_est
    phi_k_T = []; % Régresseur ligne phi(k)'
    % Termes AR (-y(k-1)...-y(k-na))
    for i = 1:na
        phi_k_T = [phi_k_T, -y_est(k-i)];
    end
    % Termes X (u(k-nk)...u(k-nb-nk+1))
    for i = 0:nb-1
         phi_k_T = [phi_k_T, u_est(k-nk-i)];
    end
    Phi_matrix_2 = [Phi_matrix_2; phi_k_T]; 
    Y_vector_2 = [Y_vector_2; y_est(k)];
end

% Calculer theta_hat_2 (LS)
try
    theta_hat_2 = (Phi_matrix_2' * Phi_matrix_2) \ (Phi_matrix_2' * Y_vector_2); 

    fprintf('\n--- Résultats Étape 10 : Identification ARX(2,2,1) par LS ---\n');
    fprintf('Paramètres estimés par LS:\n');
    fprintf('  a1 = %.6f\n', theta_hat_2(1));
    fprintf('  a2 = %.6f\n', theta_hat_2(2));
    fprintf('  b1 = %.6f\n', theta_hat_2(3));
    fprintf('  b2 = %.6f\n', theta_hat_2(4));

    % Calcul des pôles discrets (racines de A(z)=1+a1*z^-1+a2*z^-2 = 0)
    % Équivalent à chercher les racines de z^2 + a1*z + a2 = 0
    A_coeffs = [1, theta_hat_2(1), theta_hat_2(2)];
    poles_discrets = roots(A_coeffs);
    fprintf('Pôles discrets estimés: %s\n', num2str(poles_discrets.', '%.6f ')); 
    if all(abs(poles_discrets) < 1)
        fprintf('  (Le modèle est stable car les pôles sont dans le cercle unité)\n');
    else
        fprintf('  (Attention: Le modèle pourrait être instable)\n');
    end
    fprintf('------------------------------------------------------------\n');

catch ME
    fprintf('\nErreur lors du calcul LS pour ARX(2,2,1): %s\n', ME.message);
     % Arrêter ici si l'identification échoue
    return;
end

% --- Etape 11 : Comparaison des Coûts ---

% Rappel: theta_hat_2, Phi_matrix_2, Y_vector_2 calculés à l'étape 10 (pour k=3..N_est)
% Rappel: theta_hat_1 calculé à l'étape 8

% Calcul du coût pour theta_hat_2 (Modèle 2ème ordre)
Y_pred_est_2 = Phi_matrix_2 * theta_hat_2; 
Error_est_2 = Y_vector_2 - Y_pred_est_2;
Cost_LS_2 = sum(Error_est_2.^2);

fprintf('\n--- Résultats Étape 11 : Comparaison des Coûts ---\n');
fprintf('Coût J(theta_hat_2) [ARX(2,2,1)] sur données estimation (k=3..%d): %.6f\n', N_est, Cost_LS_2);

% Recalculer le coût pour theta_hat_1 (Modèle 1er ordre) sur les MÊMES points (k=3..N_est)
% Il faut reconstruire Phi pour ARX(1,1,1) mais seulement pour k=3..N_est
Phi_matrix_1_comp = []; 
Y_vector_1_comp = []; % Devrait être identique à Y_vector_2
na1=1; nb1=1; nk1=1; % Ordres pour ARX(1,1,1)
start_idx_comp = max(na, nb+nk-1)+1; % = 3, le même start_idx que pour ARX(2,2,1)

for k = start_idx_comp : N_est
    phi_k_T = []; 
    % Termes AR (-y(k-1))
    phi_k_T = [phi_k_T, -y_est(k-1)];
    % Termes X (u(k-1))
    phi_k_T = [phi_k_T, u_est(k-nk1)];
    
    Phi_matrix_1_comp = [Phi_matrix_1_comp; phi_k_T]; 
    Y_vector_1_comp = [Y_vector_1_comp; y_est(k)]; % Ce vecteur est Y_vector_2
end

% Calculer le coût de theta_hat_1 sur ces points
Y_pred_est_1_comp = Phi_matrix_1_comp * theta_hat_1; % Utilise theta_hat_1 de l'étape 8
Error_est_1_comp = Y_vector_1_comp - Y_pred_est_1_comp;
Cost_LS_1_comp = sum(Error_est_1_comp.^2);

fprintf('Coût J(theta_hat_1) [ARX(1,1,1)] sur données estimation (k=3..%d): %.6f\n', N_est, Cost_LS_1_comp);

% Comparaison
if Cost_LS_2 < Cost_LS_1_comp
    fprintf('Le modèle ARX(2,2,1) a un coût plus faible sur les données d''estimation.\n');
    fprintf('Différence de coût: %.6f\n', Cost_LS_1_comp - Cost_LS_2);
else
    fprintf('Le modèle ARX(1,1,1) a un coût plus faible ou égal (inattendu).\n');
end
fprintf('------------------------------------------------------------\n');

% --- Etape 12 : Utilisation de la fonction arx ---

% Vérifier si la toolbox est disponible (optionnel)
hasIdentToolbox = license('test', 'Identification_Toolbox');
if ~hasIdentToolbox
    warning('System Identification Toolbox non trouvée. Impossible d''exécuter cette étape.');
    return;
end

% Créer des objets iddata pour l'estimation (si pas déjà fait)
Ts = 0.1; % Rappel
data_ident_est = iddata(y_est, u_est, Ts);

% Estimer le modèle ARX(2,2,1) avec la fonction arx
na = 2; nb = 2; nk = 1;
try
    model_arx_matlab = arx(data_ident_est, [na nb nk]);

    fprintf('\n--- Résultats Étape 12 : Utilisation de arx MATLAB ---\n');
    fprintf('Modèle ARX(2,2,1) estimé par la fonction arx:\n');
    present(model_arx_matlab); % Affiche les polynômes et incertitudes

    % Extraire les coefficients pour comparaison
    A_matlab = model_arx_matlab.A; % Contient [1, a1, a2]
    B_matlab = model_arx_matlab.B; % Contient [0, b1, b2] (à cause de nk=1)
    
    % Reconstitue theta [a1, a2, b1, b2]
    theta_arx_matlab = [A_matlab(2:end), B_matlab(nk+1:nk+nb)]'; 

    fprintf('\nComparaison des paramètres (LS manuel vs arx MATLAB):\n');
    fprintf('  Param   | LS Manuel | arx MATLAB\n');
    fprintf('-------------------------------------\n');
    param_names = {'a1', 'a2', 'b1', 'b2'};
    % Rappel: theta_hat_2 contient les paramètres estimés manuellement à l'étape 10
    for i = 1:length(theta_hat_2)
        fprintf('  %-7s | %9.6f | %9.6f\n', param_names{i}, theta_hat_2(i), theta_arx_matlab(i));
    end

    % Comparaison sur les données de validation avec la fonction compare
    data_ident_val = iddata(y_val, u_val, Ts);
    figure;
    compare(data_ident_val, model_arx_matlab); 
    title('Validation du modèle ARX(2,2,1) estimé par arx MATLAB');
    grid on;
    fprintf('Le graphique "compare" montre l''ajustement (Fit) sur les données de validation.\n');
    fprintf('------------------------------------------------------------\n');

catch ME
    fprintf('\nErreur lors de l''utilisation de la fonction arx: %s\n', ME.message);
end
