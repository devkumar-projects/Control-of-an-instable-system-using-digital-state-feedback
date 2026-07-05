%% LQ control of the twin-rotor helicopter — KUMAR / DOUMBE / COUESME — Mechatronics Expertise — ENSAM
% Multivariable LQR with integral action on the 5-state helicopter model
% (elevation + pitch + yaw chain), two motor-voltage inputs, two regulated
% outputs (elevation angle and yaw rate). A spectral shift alpha_c places a
% guaranteed stability margin. The closed loop is then simulated in
% Simulink ('rebi' model) and the two outputs / two commands are plotted.

clc;
clear;
close all;

%% ==== Identified parameters ====
Kl   = -5.21894;   % yaw-channel gain
Kt   = 1.336;      % pitch-channel gain (double integrator)
Ke   = 0.29;       % elevation-channel gain
tau  = 10;         % yaw time constant (s)
xhi  = 0.08;       % elevation damping ratio
omega0 = 0.7;      % elevation natural frequency (rad/s)
n = 5;             % number of plant states
m = 2;             % number of inputs
p = 2;             % number of regulated outputs

%% ==== Plant model ====
A = [-2*xhi*omega0  -omega0^2  0     0           0;
     1                 0       0     0           0;
     0                 0       0     0           0;
     0                 0       1     0           0 ;
     0                 0       0   Kl/tau     -1/tau]

B = [Ke*omega0^2     Ke*omega0^2;
     0                     0;
     Kt                    -Kt;
     0                      0;
     0                      0]

C = [0 1 0 0 0;
     0 0 0 0 1]

% Integral augmentation: two extra states integrating the tracking errors
Aa = [A zeros(n,p); C zeros(p,p)]
Ba = [B; zeros(p,m)]

fprintf('======= System ========\n')
fprintf('---- LQ weighting ----')
Q = diag([0 0 0 0 0 40 40])   % only the integrated tracking errors are penalized
R = eye(m)                    % both motor commands weighted equally

% Spectral shift: guarantees every closed-loop pole satisfies Re(s) < -alpha_c
alphac = 0.5;
Aa_alpha = Aa + alphac * eye(n+p);

[K, P, Vp] = lqr(Aa_alpha, Ba, Q, R);

fprintf('---- Riccati solution P ----\n')
P

fprintf('---- LQ controller ----\n')
Kp = K(:, 1:n)
Ki = K(:, n+1:n+p)

fprintf('---- Closed-loop eigenvalues ----\n')
Vp

%% ==== Simulink simulation ====
fprintf('---- Running Simulink simulation (rebi model) ----\n')
sim('rebi');
t_signal = t;

%% ==== Plot: 2 outputs and 2 commands ====
figure('Name','Closed-loop numerical simulation with the LQ controller','NumberTitle','off');

subplot(4,1,1);
plot(t_signal, epsilon, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('\epsilon (rad)');
title('Elevation angle \epsilon');

subplot(4,1,2);
plot(t_signal, v, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('v (rad/s)');
title('Yaw rate v');

subplot(4,1,3);
plot(t_signal, ud, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('u_d');
title('Command u_d (right motor)');

subplot(4,1,4);
plot(t_signal, ug, 'm', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('u_g');
title('Command u_g (left motor)');

sgtitle('Closed loop with LQ correction');
