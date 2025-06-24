clc;
clear;
% === Données système continu ===
xhi = 0.0912;
w0 = 0.8038;
k = 0.3299;
T = 0.05; % Période d'échantillonnage (s)
Tc = 1; % Temps de réponse souhaité (s)
% === Modèle d'état continu ===
A = [0 1 ; -w0^2 -2*xhi*w0];
B = [0 ; k*w0];
C = [1 0];
% === Passage en discret ===
sys_disc = c2d(ss(A, B, C, 0), T, 'zoh');
[F, G, ~, ~] = ssdata(sys_disc);
% === systeme augmenté ===
Fe = [F, zeros(2,1); C, 1];
Ge = [G; 0];
% === placement de pôles) ===
Ae = [A, zeros(2,1); C, 0];
Poles_C_BO = eig(Ae);
Poles_C_BF = -1/Tc + 1i * imag(Poles_C_BO);
Poles_Z_BF = exp(Poles_C_BF * T);
% === Correcteur 1 ===
K1 = acker(Fe, Ge, Poles_Z_BF)
F_BF = Fe - Ge * K1;
B_echelon = [0; 0; -0.79]; % Consigne unitaire
C_y = [1 0 0];
C_u = -K1;
sys_y = ss(F_BF, B_echelon, C_y, 0, T);
sys_u = ss(F_BF, B_echelon, C_u, 0, T);
figure;
subplot(2,1,1);
step(sys_y, 10);
title('Réponse indicielle en sortie y en BF avec K1');
ylabel('y (rad)');

grid on;
subplot(2,1,2);
step(sys_u, 10);
title('Commande u(t)');
ylabel('u (V)');
xlabel('Temps (s)');
grid on;
figure;
hold on; grid on;
title('Placement des pôles pour Horizon à -1 et
dimensionnement de K1');
xlabel('Partie réelle'); ylabel('Partie imaginaire');
re_lim = [-3, 1];
im_lim = [-2, 2];
xlim(re_lim);
ylim(im_lim);
plot([-1/Tc -1/Tc], im_lim, 'r--', 'LineWidth', 1.2,
'DisplayName', sprintf('Horizon = %.2f', -1/Tc));
plot([0 0], im_lim, 'k', 'HandleVisibility', 'off'); % axe
imaginaire
plot(re_lim, [0 0], 'k', 'HandleVisibility', 'off'); % axe
réel
Poles_sym = -real(Poles_C_BO) + 1i * imag(Poles_C_BO);
Poles_C_BF = -1/Tc + 1i * imag(Poles_C_BO); % recalcul propre
plot(real(Poles_C_BO), imag(Poles_C_BO), 'ks',
'MarkerFaceColor', 'k', 'DisplayName', 'Pôles BO');
plot(real(Poles_sym), imag(Poles_sym), 'ko',
'MarkerFaceColor', 'none', 'DisplayName', 'Pôles Stabilisé');
plot(real(Poles_C_BF), imag(Poles_C_BF), 'o', ...
'MarkerEdgeColor', [0 0.4 0], 'MarkerFaceColor', [0 0.4
0], ...
'LineWidth', 1.3, 'DisplayName', 'Pôles en BF');
legend('Location', 'southwest');
axis equal;
hold off;

% === Ajout d'un second correcteur avec Tc2 = 0.5 s ===
Tc2 = 0.5;
Poles_C_BF2 = -1/Tc2 + 1i * imag(Poles_C_BO);
Poles_Z_BF2 = exp(Poles_C_BF2 * T);
K2 = acker(Fe, Ge, Poles_Z_BF2) % pas utilisé ici, juste pour
info
figure;
hold on; grid on;
title('Accélérations du système avec l horizon -2 et
dimensionnement de K2');
xlabel('Partie réelle'); ylabel('Partie imaginaire');
plot([0 0], [-2 2], 'k', 'HandleVisibility', 'off');
plot([-3 1], [0 0], 'k', 'HandleVisibility', 'off');
plot([-1 -1], [-2 2], 'r--', 'LineWidth', 1.2, 'DisplayName',
'Horizon Tc = 1s');
plot([-2 -2], [-2 2], 'b--', 'LineWidth', 1.2, 'DisplayName',
'Horizon Tc = 0.5s');
plot(real(Poles_C_BO), imag(Poles_C_BO), 'ks',
'MarkerFaceColor', 'k', 'DisplayName', 'Pôles BO');
plot(real(Poles_sym), imag(Poles_sym), 'ko',
'MarkerFaceColor', 'none', 'DisplayName', 'Pôles stabilisés');
plot(real(Poles_C_BF), imag(Poles_C_BF), 'o', ...
'MarkerEdgeColor', [0 0.4 0], 'MarkerFaceColor', [0 0.4
0], ...
'LineWidth', 1.3, 'DisplayName', 'Pôles BF (Tc=1s)');
plot(real(Poles_C_BF2), imag(Poles_C_BF2), 'o', ...
'MarkerEdgeColor', [0 0 0.6], 'MarkerFaceColor', [0 0 1],
...
'LineWidth', 1.3, 'DisplayName', 'Pôles BF (Tc=0.5s)');
xlim([-3 1]);
ylim([-2 2]);
axis equal;
legend('Location', 'southwest');
hold off;
F_BF2 = Fe - Ge * K2;
sys_y2 = ss(F_BF2, B_echelon, C_y, 0, T); % Système avec K2
(Tc = 0.5s)

% === réponses indicielle y2(t) ===
figure;
step(sys_y2, 10);
legend( 'Tc = 0.5s');
title('Comparaison des réponses indicielle en sortie y');
xlabel('Temps (s)');
ylabel('Élévation (rad)');
grid on;
hold off;
% === Système boucle fermée pour la commande avec le second
correcteur (Tc = 0.5s) ===
C_u2 = -K2;
sys_u2 = ss(F_BF2, B_echelon, C_u2, 0, T);
figure;
step(sys_u, 10); hold on;
step(sys_u2, 10);
legend('Commande u(t) - Tc = 1s', 'Commande u(t) - Tc =
0.5s');
title('Comparaison des commandes u(t)');
xlabel('Temps (s)');
ylabel('u (V)');
grid on;
hold off;
