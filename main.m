%{
This code implements the nonlinear output-only filtering with self-starting
integration. It estimates the states and the iputs of the system with
only response data.

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear output-only Bayesian filter for models with
              self-starting integration - Delgado Trujillo et al
%}
close all; clear all; clc
rng(1234)

%% INPUTS:
rho_infty = 1;       % high energy disipation Generalized-alpha [-] -- Page 14
tolerance = 1e-3;    % iterations tolerance                     [-] -- Page 14
max_iter  = 20;      % maximum iterations                       [-] -- Page 14

kappa = 1;    % spread parameter UT          [-] -- Page 14
alpha = 1;    % spread parameter UT          [-] -- Page 14
beta  = 1;    % non-Gaussianity parameter UT [-] -- Page 14

noise_ratio = 0.05;    % noise to signal ratio [-] -- Page 14

%{
% system
mi    = 1;      % DOFs mass                                  [1e3 kg] -- Table 1
ci    = 1;      % DOFs damping                               [kN-s/m] -- Table 1
ki    = 8;      % DOFs linear stiffness                      [kN/m]   -- Table 1
ke_dc = 3;      % DOFs EPP elastic stiffness (data creation) [kN/m]   -- Table 1
Pu_dc = 0.1;    % DOFS EPP cap strength (data creation)      [kN]     -- Table 1
ke_f  = ke_dc;  % DOFs EPP elastic stiffness (filtering)     [kN/m]   -- Page 14
Pu_f  = Pu_dc;  % DOFS EPP cap strength (filtering)          [kN]     -- Page 14
x0    = [0;     % initial condition: displacements           [m]      -- Table 1
         0;     %                    velocities              [m/s]
         0;     %                    accelerations           [m/s2]
         0];    %                    plastic deformations    [m]

% force: -- Page 14
data = load("el_centro.txt");
t    = data(:, 1);               % simulation time [s]
ft   = -9.806*mi*data(:, 2)';    % external force  [kN]

% measurements: -- Page 14
Su = 1;     % displacement selection matrix 
Sv = [];    % velocity selection matrix
Sa = 1;     % acceleration selection matrix

% filter: -- Page 14
mf_0 = 0;                % mean initial force       [kN]
Pf_0 = 1e-2;             % covariance initial force [kN2]
mx_0 = [0; 0; 0; 0];     % mean initial state       [m],  [m/s],   [m/s2],  [m]
Px_0 = (1e1)*eye(4);     % covariance initial state [m2], [m2/s2], [m2/s4], [m2]
Q    = (1e-4)*eye(4);    % covariance process noise [m2], [m2/s2], [m2/s4], [m2]
R    = [2e-3,   0;       % covar measurement noise  [m2], [m2/s4]
           0, 0.2];
E    = 0.1;              % covariance force         [kN]
%}
%{
% system
mi    = 1;      % DOFs mass                                  [1e3 kg] -- Table 1
ci    = 1;      % DOFs damping                               [kN-s/m] -- Table 1
ki    = 8;      % DOFs linear stiffness                      [kN/m]   -- Table 1
ke_dc = 3;      % DOFs EPP elastic stiffness (data creation) [kN/m]   -- Table 1
Pu_dc = 0.1;    % DOFS EPP cap strength (data creation)      [kN]     -- Table 1
ke_f  = 4;      % DOFs EPP elastic stiffness (filtering)     [kN/m]   -- Page 14
Pu_f  = 0.05;   % DOFS EPP cap strength (filtering)          [kN]     -- Page 14
x0    = [0;     % initial condition: displacements           [m]      -- Table 1
         0;     %                    velocities              [m/s]
         0;     %                    accelerations           [m/s2]
         0];    %                    plastic deformations    [m]

% force: -- Page 14
data = load("el_centro.txt");
t    = data(:, 1);               % simulation time [s]
ft   = -9.806*mi*data(:, 2)';    % external force  [kN]

% measurements: -- Page 14
Su = 1;     % displacement selection matrix
Sv = [];    % velocity selection matrix
Sa = 1;     % acceleration selection matrix

% filter: -- Page 14
mf_0 = -1;               % mean initial force       [kN]
Pf_0 = 1e-2;             % covariance initial force [kN2]
mx_0 = [0.1; 0; 0; 0.1]; % mean initial state       [m],  [m/s],   [m/s2],  [m]
Px_0 = (1e1)*eye(4);     % covariance initial state [m2], [m2/s2], [m2/s4], [m2]
Q    = (1e-4)*eye(4);    % covariance process noise [m2], [m2/s2], [m2/s4], [m2]
R    = [2e-3,   0;       % covar measurement noise  [m2], [m2/s4]
           0, 0.2];
E    = 0.1;              % covariance force         [kN]
%}
%
% system:
mi    = [1; 2; 1];         % DOFs mass                                  [1e3 kg] -- Table 2
ci    = [1; 1; 1];         % DOFs damping                               [kN-s/m] -- Table 2
ki    = [8; 4; 4];         % DOFs linear stiffness                      [kN/m]   -- Table 2
ke_dc = [3; 3; 3];         % DOFs EPP elastic stiffness (data creation) [kN/m]   -- Table 2
Pu_dc = [0.5; 0.5; 0.5];   % DOFS EPP cap strength (data creation)      [kN]     -- Table 2
ke_f  = [6; 4; 2];         % DOFs EPP elastic stiffness (filtering)     [kN/m]   -- Page 17
Pu_f  = Pu_dc;             % DOFS EPP cap strength (filtering)          [kN]     -- Page 17
x0    = zeros(12, 1);      % initial condition: displacements           [m]      -- Page 17
                           %                    velocities              [m/s]
                           %                    accelerations           [m/s2]
                           %                    plastic deformations    [m]

% force: -- Table 2
t  = (0:0.05:60)';                            % simulation time [s]
ft = [                   zeros(size(t'));     % external force  [kN]
        0.1*(2*sin(0.4*t') + sin(0.25*t'));
                    0.25*(1 - cos(0.1*t'))];

% measurements: -- Page 17
Su = [0, 1, 0];    % displacement selection matrix
Sv = [];           % velocity selection matrix
Sa = [1, 0, 0;     % acceleration selection matrix
      0, 1, 0;
      0, 0, 1];

% filter: -- Page 17
mf_0 = zeros(3, 1);      % mean initial force       [kN]
Pf_0 = 10*eye(3);        % covariance initial force [kN2]
mx_0 = 0.1*ones(12, 1);  % mean initial state       [m],  [m/s],   [m/s2],  [m]
Px_0 = 10*eye(12);       % covariance initial state [m2], [m2/s2], [m2/s4], [m2]
Q    = (1e-4)*eye(12);   % covariance process noise [m2], [m2/s2], [m2/s4], [m2]
R    = diag([0.05, ...   % covar measurement noise  [m2],
              0.5, ...   %                          [m2/s4]
                1, ...   %                          [m2/s4]
                5]);     %                          [m2/s4]
E    = diag([1e-6, ...   % covariance force         [kN]
              0.1, ...
              0.1]);
%}

%% PARAMETERS:
nt     = size(t, 1);       % number of time steps
n_Su   = size(Su, 1);      % number displacement measurements
n_Sv   = size(Sv, 1);      % number velocities measurements
n_Sa   = size(Sa, 1);      % number acceleration measurements
n_DOFs = size(mi, 1);      % number of DOFs
nx     = size(mx_0, 1);    % dimensionality of the state
nf     = size(mf_0, 1);    % dimensionality of the force

n_meas = n_Su + n_Sv + n_Sa;    % number of measurements
n_hat  = nx + nf;               % dimensionality augmented vector -- Page 11
dt     = t(2) - t(1);           % time step [s]

% weights of the UT [-] -- Eqs.15 - 17
lambda    = alpha^2*(n_hat + kappa) - n_hat;    % spread [-] -- Page 3
Wm        = NaN(2*n_hat + 1, 1);
Wm(1)     = lambda/(n_hat + lambda);            %            -- Eq.4
Wm(2:end) = 1/(2*(n_hat + lambda));             %            -- Eq.6

Wc        = NaN(2*n_hat + 1, 1);
Wc(1)     = (lambda/(n_hat + lambda)) + (1 - alpha^2 + beta);    % -- Eq.5
Wc(2:end) = 1/(2*(n_hat + lambda));                              % -- Eq.6

% measurement matrix -- Eq.14
H = zeros(n_meas, nx);
H(                1:n_Su,              1:n_DOFs) = Su;    % displacements [m]
H(       n_Su + (1:n_Sv),   n_DOFs + (1:n_DOFs)) = Sv;    % velocities    [m/s]
H(n_Su + n_Sv + (1:n_Sa), 2*n_DOFs + (1:n_DOFs)) = Sa;    % accelerations [m/s2]

% system matrices
M = sparse(diag(mi));                 % mass matrix      [1e3 kg]

C       = sparse(n_DOFs, n_DOFs);     % damping matrix   [kN - s/m]
C(1, 1) = ci(1);

K       = sparse(n_DOFs, n_DOFs);     % stiffness matrix [kN/m]
K(1, 1) = ki(1);

for i = 2:n_DOFs
    C([i - 1, i], [i - 1, i]) = C([i - 1, i], [i - 1, i]) + [ ci(i), -ci(i);
                                                             -ci(i),  ci(i)];

    K([i - 1, i], [i - 1, i]) = K([i - 1, i], [i - 1, i]) + [ ki(i), -ki(i);
                                                             -ki(i),  ki(i)];
end

theta_f  = {M; C; K; ke_f; Pu_f; rho_infty; dt;      % model parameters [-]
            tolerance; max_iter};                    % for filtering
theta_dc = {M; C; K; ke_dc; Pu_dc; rho_infty; dt;    % model parameters [-]
            tolerance; max_iter};                    % to create data

%% DATA CREATION:
% allocate space and assign the initial condition
x = NaN(nx,     nt);    x(:, 1) = x0;
y = NaN(n_meas, nt);
f = NaN(nf,     nt);

alpha_f = rho_infty/(rho_infty + 1);    % [-]

% simulate the deterministic system to get data
for i = 2:nt
    % at each time step:
    % 1) get the input of the interval -- linear approximation
    f(:, i) = (1 - alpha_f)*ft(:, i) + alpha_f*ft(:, i - 1);    % [kN]

    % 2) compute the state -- Eq.1
    x(:, i) = g(x(:, i - 1), f(:, i), theta_dc);     % [m], [m/s], [m/s2], [m]

    % 3) calculate the measurements -- Eq.13
    y(:, i) = H*x(:, i);                             % [m], [m/s], [m/s2]
end

% add noise to create the measurements
noise_mean = zeros(n_meas, 1);                            % [m], [m/s], [m/s2]
noise_sd   = sqrt(noise_ratio*var(y(:, 2:end), 0, 2));    % [m], [m/s], [m/s2]

for i = 2:nt
    y(:, i) = y(:, i) + normrnd(noise_mean, noise_sd);    % noisy measurements
end                                                       % [m], [m/s], [m/s2]

%% FILTERING:
% estimate the sytem state and the input with the nonlinear output-only filter
mx = NaN(nx,     nt);       mx(:,    1) = mx_0;
mf = NaN(nf,     nt);       mf(:,    1) = mf_0;
Px = NaN(nx, nx, nt);       Px(:, :, 1) = Px_0;
Pf = NaN(nf, nf, nt);       Pf(:, :, 1) = Pf_0;

for i = 2:nt
    % for each time step:
    % 1) predict the input -- Eqs.36 and 37
    mf_iim1 = mf(:,    i - 1);        % mean       [kN]
    Pf_iim1 = Pf(:, :, i - 1) + E;    % covariance [kN2]

    % 2) predict the state and compute the cross-covariances
    % 2.1) create the sigma points of the UT -- Eqs.50 to 56
    sqrt_Px = sqrtm(Px(:, :, i - 1));
    sqrt_Pf = sqrtm(Pf_iim1);

    Zx = NaN(nx, 2*n_hat + 1);    Zx(:, 1) = mx(:, i - 1);
    Zf = NaN(nf, 2*n_hat + 1);    Zf(:, 1) = mf_iim1;

    for j = 1:nx
        Zx(:, 1 + j)         = mx(:, i - 1) + sqrt(n_hat + lambda)*sqrt_Px(:, j);
        Zf(:, 1 + j)         = mf_iim1;
        Zx(:, 1 + j + n_hat) = mx(:, i - 1) - sqrt(n_hat + lambda)*sqrt_Px(:, j);
        Zf(:, 1 + j + n_hat) = mf_iim1;
    end

    for j = (nx + 1):n_hat
        Zx(:, 1 + j)         = mx(:, i - 1);
        Zf(:, 1 + j)         = mf_iim1 + sqrt(n_hat + lambda)*sqrt_Pf(:, j - nx);
        Zx(:, 1 + j + n_hat) = mx(:, i - 1);
        Zf(:, 1 + j + n_hat) = mf_iim1 - sqrt(n_hat + lambda)*sqrt_Pf(:, j - nx);
    end

    % 2.2) evaluate the current state (transformed sigma points) and its mean
    mx_iim1 = 0;
    g_Zhat  = NaN(nx, 2*n_hat + 1);
    for j = 1:(2*n_hat + 1)
        % current state [m], [m/s], [m/s2], [m]
        g_Zhat(:, j) = g(Zx(:, j), Zf(:, j), theta_f);

        % mean [m], [m/s], [m/s2], [m] -- Eq.37
        mx_iim1 = mx_iim1 + Wm(j)*g_Zhat(:, j);
    end

    % 2.3) compute the covariance matrixes -- Eq.39 to 41
    Px_iim1 = Q;
    Pxx     = 0;
    Pfx     = 0;
    for j = 1:(2*n_hat + 1)
        % covariance current state [m2], [m2/s2], [m2/s4], [m2]
        Px_iim1 = Px_iim1 + Wc(j)*(g_Zhat(:, j) - mx_iim1)*...
                                  (g_Zhat(:, j) - mx_iim1)';

        % cross-covariance states [m2], [m2/s2], [m2/s4], [m2]
        Pxx = Pxx + Wc(j)*(Zx(:, j) - mx(:, i - 1))*(g_Zhat(:, j) - mx_iim1)';

        % cross-covariance current state and force  [kN - m],    [kN - m/s],
        %                                           [kN - m/s2], [kN - m]
        Pfx = Pfx + Wc(j)*(Zf(:, j) - mf_iim1)*(g_Zhat(:, j) - mx_iim1)';
    end

    % 3) update the state -- Eq.42 to 46
    my          = H*mx_iim1;                      % [m],  [m/s],   [m/s2]
    Py          = H*Px_iim1*H' + R;               % [m2], [m2/s2], [m2/s4]
    Kx          = (Px_iim1*H')/Py;
    mx(:,    i) = mx_iim1 + Kx*(y(:, i)  - my);   % [m],  [m/s],   [m/s2],  [m]
    Px(:, :, i) = Px_iim1 - Kx*Py*Kx';            % [m2], [m2/s2], [m2/s4], [m2]

    % 4) update the input -- Eq.47 to 49
    Kf          = (Pfx*H')/Py;
    mf(:,    i) = mf_iim1 + Kf*(y(:, i) - my);    % [kN]
    Pf(:, :, i) = Pf_iim1 - Kf*Py*Kf';            % [kN2]
end

%% REPORT:
% get the standard deviations
x_sd = NaN(nx, nt);
f_sd = NaN(nf, nt);
for i = 1:nt
    x_sd(:, i) = sqrt(diag(Px(:, :, i)));
    f_sd(:, i) = sqrt(diag(Pf(:, :, i)));
end

for i = 1:n_DOFs
    figure('Color', 'w')
    hold on
    sgtitle(['DOF # ', num2str(i)])

    subplot(5, 1, 1)
    hold on
    patch([t; t(end:-1:1)],                       ...
          [mx(i,        :) + x_sd(i,        :),   ...
           mx(i, end:-1:1) - x_sd(i, end:-1:1)]', ...
           'r', 'FaceAlpha', 0.25, 'DisplayName', '\pm \sigma', 'EdgeColor', 'none')
    if size(Su, 1) ~= 0
        if any(Su(:, i))
            idx = find(Su(:, i));
            plot(t, y(idx, :), 'k.', 'DisplayName', 'Measurement')
        end
    end
    plot(t,  x(i, :),  'b-', 'DisplayName',      'True')
    plot(t, mx(i, :), 'r--', 'DisplayName', 'Predicted')
    axis tight
    grid on
    legend
    xlabel 't [s]'
    ylabel 'u(t) [m]'


    subplot(5, 1, 2)
    hold on
    patch([t; t(end:-1:1)],                                         ...
          [mx(n_DOFs + i,        :) + x_sd(n_DOFs + i,        :),   ...
           mx(n_DOFs + i, end:-1:1) - x_sd(n_DOFs + i, end:-1:1)]', ...
           'r', 'FaceAlpha', 0.25, 'DisplayName', '\pm \sigma', 'EdgeColor', 'none')
    if size(Sv, 1) ~= 0
        if any(Sv(:, i))
            idx = n_Su + find(Sv(:, i));
            plot(t, y(idx, :), 'k.', 'DisplayName', 'Measurement')
        end
    end
    plot(t,  x(n_DOFs + i, :),  'b-', 'DisplayName',      'True')
    plot(t, mx(n_DOFs + i, :), 'r--', 'DisplayName', 'Predicted')
    axis tight
    grid on
    legend
    xlabel 't [s]'
    ylabel 'v(t) [m/s]'


    subplot(5, 1, 3)
    hold on
    patch([t; t(end:-1:1)],                                             ...
          [mx(2*n_DOFs + i,        :) + x_sd(2*n_DOFs + i,        :),   ...
           mx(2*n_DOFs + i, end:-1:1) - x_sd(2*n_DOFs + i, end:-1:1)]', ...
           'r', 'FaceAlpha', 0.25, 'DisplayName', '\pm \sigma', 'EdgeColor', 'none')
    if size(Sa, 1) ~= 0 
        if any(Sa(:, i))
            idx = n_Su + n_Sv + find(Sa(:, i));
            plot(t, y(idx, :), 'k.', 'DisplayName', 'Measurement')
        end
    end

    plot(t,  x(2*n_DOFs + i, :),  'b-', 'DisplayName',      'True')
    plot(t, mx(2*n_DOFs + i, :), 'r--', 'DisplayName', 'Predicted')
    axis tight
    grid on
    legend
    xlabel 't [s]'
    ylabel 'a(t) [m/s2]'


    subplot(5, 1, 4)
    hold on
    patch([t; t(end:-1:1)],                                             ...
          [mx(3*n_DOFs + i,        :) + x_sd(3*n_DOFs + i,        :),   ...
           mx(3*n_DOFs + i, end:-1:1) - x_sd(3*n_DOFs + i, end:-1:1)]', ...
           'r', 'FaceAlpha', 0.25, 'DisplayName', '\pm \sigma', 'EdgeColor', 'none')
    plot(t,  x(3*n_DOFs + i, :),  'b-', 'DisplayName',      'True')
    plot(t, mx(3*n_DOFs + i, :), 'r--', 'DisplayName', 'Predicted')
    axis tight
    grid on
    legend
    xlabel 't [s]'
    ylabel '\epsilon^{(p)}(t) [m]'


    subplot(5, 1, 5)
    hold on
    patch([t; t(end:-1:1)],                       ...
          [mf(i,        :) + f_sd(i,        :),   ...
           mf(i, end:-1:1) - f_sd(i, end:-1:1)]', ...
           'r', 'FaceAlpha', 0.25, 'DisplayName', '\pm \sigma', 'EdgeColor', 'none')
    plot(t,  f(i, :),  'b-', 'DisplayName',      'True')
    plot(t, mf(i, :), 'r--', 'DisplayName', 'Predicted')
    axis tight
    grid on
    legend
    xlabel 't [s]'
    ylabel 'f(t) [kN]'
end