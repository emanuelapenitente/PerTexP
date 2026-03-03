function [S1_DFE, V1_DFE, S2_DFE, V2_DFE, R0, RC] = RC_calc(par)
%RC_calc  computes the disease-free equilibrium (DFE) and reproduction numbers (RC, R0).
%
%   INPUT
%     par : struct containing model parameters.
%
%   OUTPUT
%     S1_DFE, V1_DFE : infant susceptible / vaccinated at the DFE
%     S2_DFE, V2_DFE : non-infant susceptible / vaccinated at the DFE
%     RC            : control reproduction number (computed from NGM at the DFE)
%     R0            : basic reproduction number (computed from NGM at the DFE))

% Unpack parameters from the input struct
Lambda = par.Lambda;
r1  = par.r1; r2  = par.r2; r3  = par.r3;
eta   = par.eta; p     = par.p; psi1  = par.psi1; psi2  = par.psi2;
sigma1 = par.sigma1; sigma2 = par.sigma2;
omega = par.omega;
gamma1 = par.gamma1; gamma2 = par.gamma2;
nu = par.nu; 
beta11 = par.beta11; beta12 = par.beta12;
beta21 = par.beta21; beta22 = par.beta22;

%% Disease-Free Equilibrium (DFE) computation
% Infant block
den_S1 = 1 - r1 * (1 - psi1 - eta);
if abs(den_S1) < 1e-12
    error('Denominator for S1* is nearly zero: 1 - r1*(1 - psi1 - eta) ~ 0.');
end
S1_DFE = (1 - p) * Lambda / den_S1;

den_V1 = 1 - r1 * (1 - eta);
if abs(den_V1) < 1e-12
    error('Denominator for V1* is nearly zero: 1 - r1*(1 - eta) ~ 0.');
end
V1_DFE = (p * Lambda + r1 * psi1 * S1_DFE) / den_V1;

% Non-infant block
% Solve a 2x2 linear system for [S2*; V2*]
M = [  1 - r2 * (1 - psi2),   -r2 * omega ; ...
      -r2 * psi2,              1 - r2 * (1 - omega) ];

rhs = r1 * eta * [S1_DFE; V1_DFE];

Delta = det(M);  % = (1 - r2) * [1 - r2*(1 - psi2 - omega)]
if abs(Delta) < 1e-12
    error('Adult/non-infant system is singular: det(M) ~ 0 (check r2, psi2, omega).');
end

SV_star = M \ rhs;      % [S2*; V2*]
S2_DFE  = SV_star(1);
V2_DFE  = SV_star(2);

%% RC computation via discrete NGM: K = F * (I - T)^(-1)
% DFE population sizes
N1 = S1_DFE + V1_DFE;
N2 = S2_DFE + V2_DFE;

if N1 <= 0 || N2 <= 0
    error('DFE populations are not positive (check parameters).');
end

% Fertility (new infection) terms
c1 = r1 * (S1_DFE + sigma1 * V1_DFE);
c2 = r2 * (S2_DFE + sigma2 * V2_DFE);

% Frequency-dependent transmission scaling
g11 = beta11 / N1;
g12 = beta12 / N2;
g21 = beta21 / N1;
g22 = beta22 / N2;

Fmatrix = [ c1 * g11,  c1 * g12 ; ...
            c2 * g21,  c2 * g22 ];

% Transition matrix (infected stages progression / aging, etc.)
Tmatrix = [ r3 * (1 - gamma1 - eta),  0 ; ...
            r3 * eta,                r2 * (1 - gamma2) ];

% Check invertibility of (I - T)
A = eye(2) - Tmatrix;
if rcond(A) < 1e-12
    error('(I - T) is nearly singular: cannot compute K reliably.');
end

% Discrete Next Generation Matrix 
K  = Fmatrix / A;

% Control reproduction number: spectral radius of K
ev = eig(K);
RC = max(abs(ev));

%% R0 computation via discrete NGM: K = F * (I - T)^(-1)
F11 = r1 * beta11;
F12 = (beta12 * (1 - r2)) / eta;
F21 = (beta21 * r1 * r2 * eta) / (1 - r2);
F22 = r2 * beta22;

Fmatrix_R0 = [ F11, F12 ; ...
               F21, F22 ];

Tmatrix_R0 = [ r3 * (1 - gamma1 - eta),  0 ; ...
               r3 * eta,                r2 * (1 - gamma2) ];

A = eye(2) - Tmatrix_R0;
if rcond(A) < 1e-12
    error('(I - T) is nearly singular: cannot compute K reliably.');
end

K  = Fmatrix_R0 / A;

% Basic reproduction number: spectral radius of K
ev = eig(K);
R0 = max(abs(ev));

end
