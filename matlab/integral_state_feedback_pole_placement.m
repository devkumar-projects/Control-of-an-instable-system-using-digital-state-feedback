%% Integral state feedback by discrete pole placement — twin-rotor elevation channel
% KUMAR / DOUMBE / COUESME — Mechatronics Expertise — ENSAM
% Second-order underdamped elevation model, ZOH discretization, integral
% augmentation, Ackermann pole placement at two design horizons (Tc = 1 s
% and Tc = 0.5 s), with output and control-effort comparison.

clc;
clear;

% === Continuous-time system data (refined elevation identification) ===
xhi = 0.0912;      % damping ratio
w0  = 0.8038;      % undamped natural frequency (rad/s)
k   = 0.3299;      % static gain (rad/V)
T   = 0.05;        % sampling period (s)
Tc  = 1;           % desired closed-loop time constant (s)

% === Continuous state-space model ===
A = [0 1 ; -w0^2 -2*xhi*w0];
B = [0 ; k*w0];
C = [1 0];

% === Zero-order-hold discretization ===
sys_disc = c2d(ss(A, B, C, 0), T, 'zoh');
[F, G, ~, ~] = ssdata(sys_disc);

% === Integral-augmented system ===
Fe = [F, zeros(2,1); C, 1];
Ge = [G; 0];

% === Pole placement ===
Ae = [A, zeros(2,1); C, 0];
Poles_C_BO = eig(Ae);                          % open-loop poles (continuous)
Poles_C_BF = -1/Tc + 1i * imag(Poles_C_BO);    % desired closed-loop poles
Poles_Z_BF = exp(Poles_C_BF * T);              % mapped into the z-plane

% === Controller 1 (Tc = 1 s) ===
K1 = acker(Fe, Ge, Poles_Z_BF)
F_BF = Fe - Ge * K1;
B_step = [0; 0; -0.79];                        % unit reference injection
C_y = [1 0 0];
C_u = -K1;
sys_y = ss(F_BF, B_step, C_y, 0, T);
sys_u = ss(F_BF, B_step, C_u, 0, T);

figure;
subplot(2,1,1);
step(sys_y, 10);
title('Closed-loop step response y(t) with K1');
ylabel('y (rad)');
grid on;
subplot(2,1,2);
step(sys_u, 10);
title('Control command u(t)');
ylabel('u (V)');
xlabel('Time (s)');
grid on;

% === Pole map for the Tc = 1 s design ===
figure;
hold on; grid on;
title('Pole placement for the -1 horizon and K1 design');
xlabel('Real part'); ylabel('Imaginary part');
re_lim = [-3, 1];
im_lim = [-2, 2];
xlim(re_lim);
ylim(im_lim);
plot([-1/Tc -1/Tc], im_lim, 'r--', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Horizon = %.2f', -1/Tc));
plot([0 0], im_lim, 'k', 'HandleVisibility', 'off');   % imaginary axis
plot(re_lim, [0 0], 'k', 'HandleVisibility', 'off');   % real axis
Poles_sym  = -real(Poles_C_BO) + 1i * imag(Poles_C_BO);
Poles_C_BF = -1/Tc + 1i * imag(Poles_C_BO);            % clean recompute
plot(real(Poles_C_BO), imag(Poles_C_BO), 'ks', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Open-loop poles');
plot(real(Poles_sym), imag(Poles_sym), 'ko', ...
    'MarkerFaceColor', 'none', 'DisplayName', 'Stabilized poles');
plot(real(Poles_C_BF), imag(Poles_C_BF), 'o', ...
    'MarkerEdgeColor', [0 0.4 0], 'MarkerFaceColor', [0 0.4 0], ...
    'LineWidth', 1.3, 'DisplayName', 'Closed-loop poles');
legend('Location', 'southwest');
axis equal;
hold off;

% === Second, faster controller with Tc2 = 0.5 s ===
Tc2 = 0.5;
Poles_C_BF2 = -1/Tc2 + 1i * imag(Poles_C_BO);
Poles_Z_BF2 = exp(Poles_C_BF2 * T);
K2 = acker(Fe, Ge, Poles_Z_BF2)

figure;
hold on; grid on;
title('Accelerated design: -2 horizon and K2 sizing');
xlabel('Real part'); ylabel('Imaginary part');
plot([0 0], [-2 2], 'k', 'HandleVisibility', 'off');
plot([-3 1], [0 0], 'k', 'HandleVisibility', 'off');
plot([-1 -1], [-2 2], 'r--', 'LineWidth', 1.2, 'DisplayName', 'Horizon Tc = 1s');
plot([-2 -2], [-2 2], 'b--', 'LineWidth', 1.2, 'DisplayName', 'Horizon Tc = 0.5s');
plot(real(Poles_C_BO), imag(Poles_C_BO), 'ks', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Open-loop poles');
plot(real(Poles_sym), imag(Poles_sym), 'ko', ...
    'MarkerFaceColor', 'none', 'DisplayName', 'Stabilized poles');
plot(real(Poles_C_BF), imag(Poles_C_BF), 'o', ...
    'MarkerEdgeColor', [0 0.4 0], 'MarkerFaceColor', [0 0.4 0], ...
    'LineWidth', 1.3, 'DisplayName', 'Closed-loop poles (Tc=1s)');
plot(real(Poles_C_BF2), imag(Poles_C_BF2), 'o', ...
    'MarkerEdgeColor', [0 0 0.6], 'MarkerFaceColor', [0 0 1], ...
    'LineWidth', 1.3, 'DisplayName', 'Closed-loop poles (Tc=0.5s)');
xlim([-3 1]);
ylim([-2 2]);
axis equal;
legend('Location', 'southwest');
hold off;

F_BF2 = Fe - Ge * K2;
sys_y2 = ss(F_BF2, B_step, C_y, 0, T);   % system with K2 (Tc = 0.5 s)

% === Step response y2(t) ===
figure;
step(sys_y2, 10);
legend('Tc = 0.5s');
title('Output step-response comparison');
xlabel('Time (s)');
ylabel('Elevation (rad)');
grid on;
hold off;

% === Closed-loop control effort with the second controller (Tc = 0.5 s) ===
C_u2 = -K2;
sys_u2 = ss(F_BF2, B_step, C_u2, 0, T);
figure;
step(sys_u, 10); hold on;
step(sys_u2, 10);
legend('Command u(t) - Tc = 1s', 'Command u(t) - Tc = 0.5s');
title('Control-effort comparison u(t)');
xlabel('Time (s)');
ylabel('u (V)');
grid on;
hold off;
