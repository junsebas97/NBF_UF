function [mx, Px, mfu, Pfu] = this_filter(y, mx_0, Px_0, mfu_0, Pfu_0, E, Q, ...
                                          R, kappa, alpha, beta, theta)
%{
this function is the proposed nonlinear filter with uncertain forces. It
infers the states and the uncertain forces of the given system using ramdon
walks, the Unscented Transfrom (UT), and SSTIAs.

INPUTS:
y:     measurements                              [m],  [m/s],   [m/s2]
mx_0:  initial state mean                        [m],  [m/s],   [m/s2],  [m]
Px_0:  initial state covariance                  [m2], [m2/s2], [m2/s4], [m2]
mfu_0: mean of the initial uncertain force       [kN]
Pfu_0: covariance of the initial uncertain force [kN2]
E:     covariance of uncertain force random walk [kN2]
Q:     covariance process noise                  [m2], [m2/s2], [m2/s4], [m2]
R:     covariance measurement noise              [m2], [m2/s2], [m2/s4]
kappa: spread parameter UT                       [-]
alpha: spread parameter UT                       [-]
beta:  non-Gaussianity parameter UT              [-]
theta: system parameters. They must organized as
    ---                                       ---
   | measurement matrix                 [-]      |
   | known forces                       [kN]     |
   | influence vector unknown forces    [-]      |
   | influence vector uncertain forces  [-]      |
   | mass matrix                        [1e3 kg] |
   | damping matrix                     [kN-s/m] |
   | stiffnesses                        [kN/m]   |
   | post-yield stiffness ratio         [-]      |
   | parameter of the Bouc-Wen model    [-]      |
   | parameter of the Bouc-Wen model    [-]      |
   | parameter of the Bouc-Wen model    [-]      |
   | high energy dissipation            [-]      |
   | time step                          [s]      |
   | tolerance Newton-Raphson           [-]      |
   | maximum iterations Newton-Raphson  [-]      |
   ---                                        ---

OUTPUTs:
mx:  state mean                     [m],  [m/s],   [m/s2],  [m]
Px:  state covariance               [m2], [m2/s2], [m2/s4], [m2]
mfu: mean of uncertain forces       [kN]
Pfu: covariance of uncertain forces [kN]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}
%% PARAMETRIZATION:
% extract the system parameters
H         = theta{1};          % measurement matrix                    [-]
f_k       = theta{2};          % known forces                          [kN]
S_fk      = theta{3};          % influence vector for unknown forces   [-]
S_fu      = theta{4};          % influence vector for uncertain forces [-]
rho_infty = theta{12};         % high energy dissipation of the TIA    [-]
theta     = theta(5:15);       % system parameters for process model

Nt        = size(f_k, 2);      % number of time steps
nx        = size(mx_0, 1);     % dimensionality of the state
nfu       = size(mfu_0, 1);    % dimensionality of the uncerain force
nz        = nx + nfu;          % augmented dimensionality -- {1} Page 4

% evaluate the known force at the sampling time
alpha_f    = rho_infty/(rho_infty + 1);          % force factor [-] -- {2} Eq.25
f_k(2:end) = alpha_f*f_k(:, (2:end) - 1) + ...   % known force  [kN]
             (1 - alpha_f)*f_k(:, (2:end));

% compute the weights of the UT [-]
lambda    = alpha^2*(nz + kappa) - nz;                       % [-] -- {1} Page 3
Wm        = NaN(2*nz + 1, 1);
Wc        = NaN(1, 2*nz + 1);
Wm(1)     = lambda/(nz + lambda);                            % [-] -- {1} Eq.7
Wc(1)     = (lambda/(nz + lambda)) + (1 - alpha^2 + beta);   % [-] -- {1} Eq.8
Wm(2:end) = 1/(2*(nz + lambda));                             % [-] -- {1} Eq.9
Wc(2:end) = 1/(2*(nz + lambda));                             % [-] -- {1} Eq.9

%% FILTERING:
% apply the filter to estimate the sytem state and the uncertain force
mx  = NaN( nx,      Nt);       mx( :,    1) = mx_0;
mfu = NaN(nfu,      Nt);       mfu(:,    1) = mfu_0;
Px  = NaN( nx,  nx, Nt);       Px( :, :, 1) = Px_0;
Pfu = NaN(nfu, nfu, Nt);       Pfu(:, :, 1) = Pfu_0;

for i = 2:Nt
    % for each time step:
    % 1) predict the input
    mfu_iim1 = mfu(:,    i - 1);        % mean       [kN]  -- {1} Eq.34
    Pfu_iim1 = Pfu(:, :, i - 1) + E;    % covariance [kN2] -- {1} Eq.35

    % 2) predict the state and compute its covariances
    % 2.1) create the sigma points of the UT
    mz      = [mx(:, i - 1); mfu_iim1];                         %      {1} Eq.20
    sqrt_Pz = sqrt(nz + lambda)*...                             % from {1} Eq.21
              chol(blkdiag(Px(:, :, i - 1), Pfu_iim1), 'lower');

    Z                     = NaN(nz, 2*nz + 1);
    Z(:,               1) = mz;                   % {1} Eq.1
    Z(:, 1 +      (1:nz)) = mz + sqrt_Pz;         % {1} Eq.2
    Z(:, 1 + nz + (1:nz)) = mz - sqrt_Pz;         % {1} Eq.3

    X = Z(      (1:nx), :);         % [m], [m/s], [m/s2], [m] -- {1} Page 4
    U = Z(nx + (1:nfu), :);         % [kN]                    -- {1} Page 4
    F = S_fk*f_k(:, i) + S_fu*U;    % [kN]                    -- {1} Page 4

    % 2.2) evaluate the current state by propagating the sigma points
    g_XF = NaN(nx, 2*nz + 1);
    for j = 1:(2*nz + 1)
        g_XF(:, j) = g(X(:, j), F(:, j), theta);   % [m], [m/s], [m/s2], [m]
    end

    % 2.3) compute the statistics of the current state
    mx_iim1 = g_XF*Wm;                        % mean       [m],     -- {1} Eq.36
                                              %            [m/s], 
                                              %            [m/s2],
                                              %            [m]
    g_dist  = g_XF - mx_iim1;
    Px_iim1 = (Wc.*g_dist)*g_dist' + Q;       % covariance [m2],    -- {1} Eq.37
                                              %            [m2/s2],
                                              %            [m2/s4],
                                              %            [m2]
    
    % 3) update the state and uncertain force -- {1} Eqs.44 to 52
    my           = H*mx_iim1;                     % [m],    [m/s],    [m/s2]
    Py           = H*Px_iim1*H' + R;              % [m2],   [m2/s2],  [m2/s4]
    P_fux        = (Wc.*(U - mfu_iim1))*g_dist';  % [kN-m], [kN-m/s], [kN-m/s2], [kN-m]

    Kf           = (P_fux*H')/Py;
    mfu(:,    i) = mfu_iim1 + Kf*(y(:, i) - my);  % [kN]
    Pf_aux       = Pfu_iim1 - Kf*Py*Kf';          % [kN2]
    Pfu(:, :, i) = (1/2)*(Pf_aux + Pf_aux');

    Kx           = (Px_iim1*H')/Py;
    mx(:,    i)  = mx_iim1 + Kx*(y(:, i)  - my);  % [m],  [m/s],   [m/s2],  [m]
    Px_aux       = Px_iim1 - Kx*Py*Kx';           % [m2], [m2/s2], [m2/s4], [m2]
    Px(:, :, i)  = (1/2)*(Px_aux + Px_aux');
end
end

function [x_i] = g(x_im1, f_i, theta)
%{
this function is the SSTIA-built process model. It computes the state of the
given system with the previous state and the current force using the
generalised-alpha method. The system uses the Bouc-Wen model for the
restoring force.

INPUTS:
x_im1: previous system state  [m], [m/s], [m/s2], [m]
f_i:   current force          [kN]
theta: system parameters.

OUTPUTs:
x_i: current system state  [m], [m/s], [m/s2], [m]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}

% get the system properties
M         = theta{1};      % mass matrix                       [1e3 kg]
C         = theta{2};      % damping matrix                    [kN-s/m]
k         = theta{3};      % stiffnesses                       [kN/m]
alpha_BW  = theta{4};      % post-yield stiffness ratio        [-]
beta_BW   = theta{5};      % parameter of the Bouc-Wen model   [-]
gamma_BW  = theta{6};      % parameter of the Bouc-Wen model   [-]
n_BW      = theta{7};      % parameter of the Bouc-Wen model   [-]
rho_infty = theta{8};      % high energy dissipation           [-]
dt        = theta{9};      % time step                         [s]
tolerance = theta{10};     % tolerance Newton-Raphson          [-]
max_iter  = theta{11};     % maximum iterations Newton-Raphson [-]

N_DOFs    = size(M, 1);    % numbers of DOFs

% calculate the parameters of the generalised-alpha method
alpha_f = rho_infty/(rho_infty + 1);            % [-] -- {2} Eq.25
alpha_m = (2*rho_infty - 1)/(rho_infty + 1);    % [-] -- {2} Eq.25
gamma   = 1/2 - alpha_m + alpha_f;              % [-] -- {2} Eq.17
beta    = (1/4)*(1 - alpha_m + alpha_f)^2;      % [-] -- {2} Eq.20

% extract the previous state -- {1} Page 3
u_tim1  = x_im1(           (1:N_DOFs));    % displacement      [m]
v_tim1  = x_im1(  N_DOFs + (1:N_DOFs));    % velocity          [m/s]
a_tim1  = x_im1(2*N_DOFs + (1:N_DOFs));    % acceleration      [m/s2]
xi_tim1 = x_im1(3*N_DOFs + (1:N_DOFs));    % BW displacements  [m]

% assess the previous restoring force
r_tim1 = alpha_BW.*k.*diff([0; u_tim1]) + ...    % springs force [kN] 
         (1 - alpha_BW).*k.*xi_tim1;             % -- {1} Eq.54
r_tim1 = r_tim1 - [r_tim1(2:end); 0];            % DOFs force    [kN]

% initalize the displacements and velocities
u_ti = u_tim1;    % [m]
v_ti = v_tim1;    % [m/s]

% compute the equivalent stiffness matrix and force vector
K_hat = M*((1 - alpha_m)/(beta*dt^2))     +                      ... % [kN/m]
        C*((1 - alpha_f)*gamma/(beta*dt));

f_hat = f_i                                                    - ... % [kN]
        (  M*(1 - 1/(2*beta) + alpha_m/(2*beta))               + ...
           C*(dt*(1 - alpha_f)*(1 - gamma/(2*beta)))  )*a_tim1 - ...
        (  C*(1 - gamma/beta + alpha_f*gamma/beta)             - ...
           M*((1 - alpha_m)/(beta*dt))                )*v_tim1 - ...
        (- M*((1 - alpha_m)/(beta*dt^2))                       - ...
           C*((1 - alpha_f)*gamma/(beta*dt))          )*u_tim1 - ...
        alpha_f*r_tim1;

% find the new nodal displacements with Newton-Raphson iteration
for iter = 1:max_iter
    % in each iteration,
    % 1) calculate the BW response
    [xi_ti, r_ti, dr_du] = BW_model(xi_tim1, v_tim1, v_ti, u_ti, ...   % [m],
                                    k, alpha_BW, beta_BW,        ...   % [kN],
                                    gamma_BW, n_BW, dt);               % [kN/m],

    % 2) assess the tangent stiffness and compute the residual force
    K_t   = K_hat + (1 - alpha_f)*dr_du;                  % [kN/m]
    R_for = f_hat - (K_hat*u_ti + (1 - alpha_f)*r_ti);    % [kN]

    % 3) if tolerance is not met, update the current displacements and
    %    velocities
    if norm(R_for) > tolerance
        duN  = K_t\R_for;                                 % increment    [m]
        u_ti = u_ti + duN;                                % displacement [m]
        vt_i = (gamma/(beta*dt))*(u_ti - u_tim1) + ...    % velocity     [m/s]
               (1 - gamma/beta)*v_tim1           + ... 
               dt*(1 - gamma/(2*beta))*a_tim1;
    else
        break
    end
end

% compute the current accelerations and form the current state vector
at_i = (1/(beta*dt^2))*(u_ti - u_tim1) - ...    % acceleration [m/s2] 
       (1/(beta*dt))*v_tim1            - ...
       (1/(2*beta) - 1)*a_tim1;

x_i  = [u_ti; vt_i; at_i; xi_ti];               % state [m], [m/s], [m/s2], [m]
end                                             %  -- {1} Page 3