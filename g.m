function [x_i] = g(x_im1, f_i, theta)
%{
this function is the self-starting TIA state function. It computes the state of
the given system with the previous state and the current input using the
Generalized-alpha method. The system uses an Elastic-Perfectly Plastic (EPP)
model for the nonlinear part of the restoring force.

INPUTS:
xt_im1: previous system state:
    ---                            ---
   | displacements            [m]    |
   | acceleration             [m/s2] |
   | velocities               [m/s]  |
   | EEP plastic deformations [m]    |
   ---                             ---

ft_i:   current step external force [kN]
theta: system parameters:
    ---                                      ---
   | mass matrix                       [1e3 kg] |
   | damping matrix                    [kN-s/m] |
   | elastic stiffness matrix          [kN/m]   |
   | initial stiffness EPP             [kN/m]   |
   | cap strengths EPP                 [kN]     |
   | high energy dissipation           [-]      |
   | time step                         [s]      |
   | tolerance Newton-Raphson          [-]      |
   | maximum iterations Newton-Raphson [-]      |
   ---                                       ---

OUTPUTs:
xt_i: current system state:
    ---                            ---
   | displacements            [m]    |
   | acceleration             [m/s2] |
   | velocities               [m/s]  |
   | EEP plastic deformations [m]    |
   ---                             ---

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} An output-only Bayesian filter for nonlinear models with
                  Generalized-alpha integration method - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}

%% PARAMETRIZATION:
% get the system properties
M         = theta{1};      % mass matrix                       [1e3 kg]
C         = theta{2};      % damping matrix                    [kN-s/m]
K         = theta{3};      % elastic stiffness matrix          [kN/m]
ki        = theta{4};      % initial stiffness EPP             [kN/m]
Pu        = theta{5};      % cap strengths EPP                 [kN]
rho_infty = theta{6};      % high energy dissipation           [-]
dt        = theta{7};      % time step                         [s]
tolerance = theta{8};      % tolerance Newton-Raphson          [-]
max_iter  = theta{9};      % maximum iterations Newton-Raphson [-]

n_DOFs    = size(M, 1);    % numbers of DOFs

% calculate the parameters of the TIA
alpha_f = rho_infty/(rho_infty + 1);            % [-] -- {2} Eq.25
alpha_m = (2*rho_infty - 1)/(rho_infty + 1);    % [-] -- {2} Eq.25
gamma   = 1/2 - alpha_m + alpha_f;              % [-] -- {2} Eq.17
beta    = (1/4)*(1 - alpha_m + alpha_f)^2;      % [-] -- {2} Eq.20

% extract the previous state -- {1} Page 2
u_tim1  = x_im1(             1:  n_DOFs);    % displacement         [m]
v_tim1  = x_im1(  (n_DOFs + 1):2*n_DOFs);    % velocity             [m/s]
a_tim1  = x_im1((2*n_DOFs + 1):3*n_DOFs);    % acceleration         [m/s2]
xi_tim1 = x_im1((3*n_DOFs + 1):     end);    % plastic deformations [m]

%% PREVIOUS RESPONSE:
% assess the previous force of the EPP elements
[p_tim1, ~, ~] = EPP_model(u_tim1, xi_tim1, Pu, ki);    % [kN]

%% TIME-INTEGRATION:
% initalize the displacements
u_ti  = u_tim1;     % DOFs displacements    [m]
xi_ti = xi_tim1;    % plastic displacements [m]

% compute the equivalent stiffness matrix and force vector
K_hat = M*((1 - alpha_m)/(beta*dt^2))     +                      ... % [kN/m]
        C*((1 - alpha_f)*gamma/(beta*dt)) + K*(1 - alpha_f);

f_hat = f_i                                                   - ... % [kN]
        (  M*(1 - 1/(2*beta) + alpha_m/(2*beta))               + ...
           C*(dt*(1 - alpha_f)*(1 - gamma/(2*beta)))  )*a_tim1 - ...
        (  C*(1 - gamma/beta + alpha_f*gamma/beta)             - ...
           M*((1 - alpha_m)/(beta*dt))                )*v_tim1 - ...
        (  alpha_f*K - M*((1 - alpha_m)/(beta*dt^2))           - ...
           C*((1 - alpha_f)*gamma/(beta*dt))          )*u_tim1 - ...
        alpha_f*p_tim1;

% find the new nodal displacements with Newton-Raphson iteration
for iter = 1:max_iter
    % 1) calculate the EPP response
    [pi, dpi_du, xi_ti] = EPP_model(u_ti, xi_ti, Pu, ki);    % [kN], [kN/m], [m]

    % 2) assess the tangent stiffness and compute the residual force
    K_t   = K_hat + (1 - alpha_f)*dpi_du;               % [kN/m]
    R_for = f_hat - (K_hat*u_ti + (1 - alpha_f)*pi);    % [kN]

    % 3) if tolerance is not met, update the current displacements
    if norm(R_for) > tolerance
        duN  = K_t\R_for;      % increment [m]
        u_ti = u_ti + duN;     %           [m]
    else
        break
    end
end

% 4) compute the current velocities and accelerations
vt_i = (gamma/(beta*dt))*(u_ti - u_tim1) +                       ...    % [m/s]
        (1 - gamma/beta)*v_tim1 + dt*(1 - gamma/(2*beta))*a_tim1;

at_i = (1/(beta*dt^2))*(u_ti - u_tim1) - (1/(beta*dt))*v_tim1 -  ...    % [m/s2]
       (1/(2*beta) - 1)*a_tim1;

%% FORMAT:
% form the current state vector -- {1} Page 2
x_i = [u_ti; vt_i; at_i; xi_ti];    % [m], [m/s], [m/s2], [m]

end