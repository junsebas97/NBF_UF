function [p, Kn, eps_p] = EPP_model(u, eps_p, Pu, ke)
%{
this function is the EPP force model. it returns the restoring force vector and
tangent stiffness matrix.

INPUTS:
u:  displacements                            [m]
up: plastic displacement of the EPP elements [m]
Pu: cap strength of the EPP elements         [kN]
ke: initial stiffness of the EPP elements    [kN/m]

OUTPUTS:
p:  EPP restoring force vector                       [kN]
Kn: EPP tangent stiffness matrix                     [kN/m]
up: updated plastic displacement of the EPP elements [m]

MADE BY:      junsebas97
BIBLIOGRAPHY: An output-only Bayesian filter for nonlinear models with
              Generalized-alpha integration method - Delgado Trujillo et al
%}
n_EPP = size(Pu, 1);               % number of DOFs
p     = sparse(n_EPP, 1);
Kn    = sparse(n_EPP, n_EPP);

for i = 1:n_EPP
    % in each EPP element:
    % 1) obtain the deformation and estimate the force,
    if i == 1;    eps = u(i);               % [m]
    else;         eps = u(i) - u(i - 1);    % [m]  -- Page 12
    end

    pi = ke(i)*(eps - eps_p(i));            % [kN] -- Page 12

    % 2) compute the nonlinear stiffness, and correct the force and plastic
    %    displacements if it does not remain elastic,
    if norm(pi) <= Pu(i)
        dpdu = ke(i);    % nonlinear stiffness [kN/m] -- Eq.57
    else
        ey_i     = Pu(i)/ke(i);                        % [m] -- Page 14
        de_p     = (eps - eps_p(i))*...                % [m] -- Eq.59
                   (1 - ey_i/norm(eps - eps_p(i)));      
        eps_p(i) = eps_p(i) + de_p;                    % [m] -- Eq.58

        dpdu = 0;                         % nonlinear stiffness [kN/m] -- Eq.57
        pi   = ke(i)*(eps - eps_p(i));    % corrected force     [kN]   -- Eq.60
    end

    % 3) translate the response to the force vector and the stiffness
    % matrix -- Eqs.61 and 62
    if i == 1
        p(i,  1) = p(i, 1)  + pi;      % [kN]
        Kn(i, i) = Kn(i, i) + dpdu;    % [kN/m]
    else
        p( (i - 1):i,         1) = p((i - 1):i, 1) + [pi; -pi];      % [kN]
        Kn((i - 1):i, (i - 1):i) = Kn((i - 1):i, (i - 1):i) + ...    % [kN/m]
                                   [ dpdu, -dpdu;
                                    -dpdu,  dpdu];
    end
end
end