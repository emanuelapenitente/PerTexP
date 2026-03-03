function [x, G1, G2, inc1, inc2, vax1, vax2, deaths, R1, R2] = discrete_system(t0, t1, x0, par)
% Discrete-time pertussis simulator with 6 state variables:
%         (S1, V1, I1, S2, V2, I2), with a weekly time step.
%
% The system is iterated from time step t0 to time step t1 (included).
%
% INPUT
%   t0, t1 : integers with t1 > t0
%   x0     : [6x1] initial conditions [S1; V1; I1; S2; V2; I2]
%   par    : struct containing model parameters
%
% OUTPUT
%   x      : [6 x (t1-t0+1)] solution trajectory (column 1 = x0)
%   G1,G2  : [1 x (t1-t0)] infection probabilities over each step
%   inc1,inc2 : [1 x (t1-t0)] incidences (new infections per step)
%   vax1,vax2 : [1 x (t1-t0)] vaccinations (new doses per step)
%   deaths : [1 x (t1-t0)] disease-induced deaths per step (infants)
%   R1,R2  : [1 x (t1-t0+1)] recovered compartments (computed by mass balance)

%Allocation and initialization
T = t1 - t0;        % number of discrete time steps
x = zeros(6, T + 1);
x(:, 1) = x0(:);

G1     = zeros(1, T);      G2     = zeros(1, T);
inc1   = zeros(1, T);      inc2   = zeros(1, T);
vax1   = zeros(1, T);     vax2   = zeros(1, T);
deaths = zeros(1, T);

R1 = zeros(1, T + 1);    R2 = zeros(1, T + 1);

% Unpack parameters from struct
N1 = par.N1; N2 = par.N2;
Lambda = par.Lambda; 
r1 = par.r1; r2 = par.r2; r3 = par.r3; d = par.d; 
eta = par.eta;
p = par.p; 
psi1 = par.psi1; psi2 = par.psi2;
sigma1 = par.sigma1; sigma2 = par.sigma2;
omega = par.omega; 
gamma1 = par.gamma1; gamma2 = par.gamma2; 
nu = par.nu;
beta11 = par.beta11; beta12 = par.beta12;
beta21 = par.beta21; beta22 = par.beta22;

% Simulation
for t = 1:T
    % Current state
    S1 = x(1, t);  V1 = x(2, t);  I1 = x(3, t);
    S2 = x(4, t);  V2 = x(5, t);  I2 = x(6, t);

    % Recovered compartments (not state variables, computed by mass balance)
    R1(t) = N1 - S1 - V1 - I1;
    R2(t) = N2 - S2 - V2 - I2;

    % Infection probabilities over the interval [t, t+1] (i.e., probability of becoming infected during one time step)
    G1(t) = 1 - exp(-(beta11 * (I1 / N1) + beta12 * (I2 / N2)));
    G2(t) = 1 - exp(-(beta21 * (I1 / N1) + beta22 * (I2 / N2)));

    % Incidences (new infections from t to t+1)
    inc1(t) = r1 * G1(t) * (S1 + sigma1 * V1);
    inc2(t) = r2 * G2(t) * (S2 + sigma2 * V2);

    % Vaccinations (new doses from t to t+1)
    vax1(t) = p * Lambda + r1 * psi1 * S1;
    vax2(t) = r2 * psi2 * S2;

    % Disease-induced deaths (new deaths from t to t+1)
    deaths(t) = d * I1;

    % State update at the next time step
    % Group 1 (infants)
    S1_next = (1 - p) * Lambda  + r1 * (1 - G1(t)) * S1 - r1 * psi1 * S1 - r1 * eta * S1;
    V1_next = p * Lambda  + r1 * (1 - sigma1 * G1(t)) * V1  + r1 * psi1 * S1  - r1 * eta * V1;
    I1_next = r1 * G1(t) * (S1 + sigma1 * V1) + r3 * (1 - gamma1) * I1  - r3 * eta * I1;

    % Group 2 (non-infants)
    S2_next = r1 * eta * S1 + r2 * nu * R2(t) + r2 * (1 - G2(t)) * S2 - r2 * psi2 * S2 + r2 * omega * V2;
    V2_next = r1 * eta * V1 + r2 * psi2 * S2 - r2 * omega * V2 + r2 * (1 - sigma2 * G2(t)) * V2;
    I2_next = r3 * eta * I1 + r2 * G2(t) * (S2 + sigma2 * V2) + r2 * (1 - gamma2) * I2;

    % Store updated state
    x(:, t + 1) = [S1_next; V1_next; I1_next; S2_next; V2_next; I2_next];
end

% Final recovered compartments (at t = T+1)
R1(T + 1) = N1 - x(1, T + 1) - x(2, T + 1) - x(3, T + 1);
R2(T + 1) = N2 - x(4, T + 1) - x(5, T + 1) - x(6, T + 1);

end

