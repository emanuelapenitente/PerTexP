function [S1_DFE, V1_DFE, S2_DFE, V2_DFE, R0, RC] = RC_calc(par)
%RC_CALC Compute the disease-free equilibrium and reproduction numbers.

% INPUT
%   par : model parameter structure
%
% OUTPUT
%   S1_DFE, V1_DFE : infant susceptible/vaccinated states at the DFE
%   S2_DFE, V2_DFE : non-infant susceptible/vaccinated states at the DFE
%   R0              : basic reproduction number (no vaccination)
%   RC              : control reproduction number at the vaccinated DFE

% Unpack parameters
N1 = par.N1;
N2 = par.N2;

Lambda = par.Lambda;
Mnet = par.M;

r1 = par.r1;
r2 = par.r2;
r3 = par.r3;
eta = par.eta;

p = par.p;
psi1 = par.psi1;
psi2 = par.psi2;
sigma1 = par.sigma1;
sigma2 = par.sigma2;
omega = par.omega;

gamma1 = par.gamma1;
gamma2 = par.gamma2;

beta11 = par.beta11;
beta12 = par.beta12;
beta21 = par.beta21;
beta22 = par.beta22;

%% Disease-free equilibrium: infant block
den_S1 = 1 - r1 * (1 - psi1 - eta);
if abs(den_S1) < 1e-12
    error('Denominator for S1* is nearly zero.');
end
S1_DFE = (1 - p) * Lambda / den_S1;

den_V1 = 1 - r1 * (1 - eta);
if abs(den_V1) < 1e-12
    error('Denominator for V1* is nearly zero.');
end
V1_DFE = (p * Lambda + r1 * psi1 * S1_DFE) / den_V1;

%% Disease-free equilibrium: non-infant block
A_SV = [ 1 - r2 * (1 - psi2),   -r2 * omega ; ...
        -r2 * psi2,               1 - r2 * (1 - omega) ];

rhs = [Mnet + r1 * eta * S1_DFE; ...
                r1 * eta * V1_DFE];

if rcond(A_SV) < 1e-12
    error('Non-infant DFE linear system is nearly singular.');
end

SV_star = A_SV \ rhs;
S2_DFE = SV_star(1);
V2_DFE = SV_star(2);

% Demographic consistency checks.  The model forces of infection use the
% prescribed population sizes par.N1 and par.N2, so the DFE should reproduce
% those totals up to numerical precision.
if abs((S1_DFE + V1_DFE) - N1) > 1e-6 * max(1,N1)
    warning('Infant DFE total differs from N1 by %.6g.', ...
        (S1_DFE + V1_DFE) - N1);
end
if abs((S2_DFE + V2_DFE) - N2) > 1e-6 * max(1,N2)
    warning('Non-infant DFE total differs from N2 by %.6g.', ...
        (S2_DFE + V2_DFE) - N2);
end

%% Control reproduction number RC
% Discrete next-generation matrix K = F*(I-T)^(-1).

c1 = r1 * (S1_DFE + sigma1 * V1_DFE);
c2 = r2 * (S2_DFE + sigma2 * V2_DFE);

Fmatrix = [ c1 * beta11 / N1,  c1 * beta12 / N2 ; ...
            c2 * beta21 / N1,  c2 * beta22 / N2 ];

Tmatrix = [ r3 * (1 - gamma1 - eta),  0 ; ...
            r3 * eta,                 r2 * (1 - gamma2) ];

A = eye(2) - Tmatrix;
if rcond(A) < 1e-12
    error('(I-T) is nearly singular: cannot compute RC reliably.');
end

K = Fmatrix / A;
RC = max(abs(eig(K)));

%% Basic reproduction number R0


Fmatrix_R0 = [ r1 * beta11,                    r1 * (N1/N2) * beta12 ; ...
               r2 * (N2/N1) * beta21,          r2 * beta22 ];

A0 = eye(2) - Tmatrix;
if rcond(A0) < 1e-12
    error('(I-T) is nearly singular: cannot compute R0 reliably.');
end

K0 = Fmatrix_R0 / A0;
R0 = max(abs(eig(K0)));

end
