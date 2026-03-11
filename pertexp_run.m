function out = pertexp_run(gui)
%PERTEXP_RUN  Core PerTexP simulator (called by the GUI).
%
% INPUT
%   gui : struct with the information given by the app. Expected fields:
%         - gui.p           : maternal vaccination coverage (fraction in [0,1])
%         - gui.psi1_ann    : annual infant DTaP vaccination probability (fraction in [0,1])
%         - gui.psi2_ann    : annual non-infant Tdap booster probability (fraction in [0,1])
%         - gui.sigma1      : vaccine failure in infants  (0 = perfect protection)
%         - gui.sigma2      : vaccine failure in non-infants
%         - gui.vax_imm_year: duration (years) of vaccine-induced immunity (used also for natural immunity here)
%         - gui.x0          : initial state vector [S1; V1; I1; S2; V2; I2]
%         - gui.noiseOn     : logical flag enabling/disabling inter-annual state noise
%
% OUTPUT
%   out : struct with the quantities needed by the GUI barplots:
%         - out.RC, out.R0
%         - out.inc1_annual, out.inc2_annual
%         - out.vax1_annual, out.vax2_annual
%         - out.year_vector

%% --- Model parameters ---
% --- Demography / population sizes (Italy, January 2024) ---
N  = 58971230;              % total population
N1 = 380630;                % infants (< 1 year)
N2 = N - N1;                % non-infants (>= 1 year)

% --- Ageing from infants to non-infants (weekly probability) ---
eta_ann = 1/52;             % annual rate (1 year^-1)
eta     = 1 - exp(-eta_ann);% weekly probability of aging out of infants

% --- Survival probabilities, infants ---
mu1_ann_per_thousand = 2.57230;
mu1_ann              = mu1_ann_per_thousand/1000;
r1_ann               = 1 - mu1_ann_per_thousand/10000;  
r1                   = r1_ann^(1/52);

% --- Survival probabilities, infants ---
mu1_ann_per_thousand = 2.57230;                   % annual death probability per 1000 live births
mu1_ann              = mu1_ann_per_thousand/1000; % annual death probability
r1_ann               = 1 - mu1_ann_per_thousand/10000; % annual survival probability
r1                   = r1_ann^(1/52);             % weekly survival probability

% --- Calibration to match desired equilibria ---
Lambda = (1 - r1*(1 - eta)) * N1;                 % weekly births into infants
r2     = 1 - (r1 * eta * N1) / N2;                % weekly survival probability, non-infants
mu2_ann = 1 - r2^52;                              % annual death probability, non-infants (derived)
nu_rate = 1/(10*52);                              % waning rate of natural immunity (weeks^-1)
nu      = 1 - exp(-nu_rate);                      % weekly probability of losing natural immunity

% --- Recovery (weekly probabilities) ---
gamma1_rate = 1/3;                % recovery rate (infants)
gamma1      = 1 - exp(-gamma1_rate); % weekly recovery probability (infants)

gamma2_rate = 1/3;                % recovery rate (non-infants)
gamma2      = 1 - exp(-gamma2_rate); % weekly recovery probability (non-infants)

% --- Disease-induced mortality, infants ---
p_CFR  = 6.3/1000;                                %Case Fatality Rate (CFR)
r3 = (1 - p_CFR)/((1 - p_CFR) + p_CFR*(gamma1 + eta)); % reduced survival rate for infected infants
if r3 <= 0 || r3 >= 1
    warning('Computed r3 is outside (0,1): r3 = %.6f', r3);
end
d = r1-r3;          % weekly disease-induced death probability


% --- Transmission rates ---
beta11 = 0.81;    beta12 = 3.0;
beta21 = 0.002;   beta22 = 0.3;

%% --- Vaccination parameters (GUI inputs + unit conversions) ---
% The GUI provides annual quantities and/or user-friendly units.
% Here we convert everything to weekly probabilities as required by the simulator.
p            = gui.p;             % already a fraction (0-1)
psi1_ann     = gui.psi1_ann;      % annual probability (fraction in [0,1])
psi2_ann     = gui.psi2_ann;      % annual probability (fraction in [0,1], already converted from "per 10,000")
sigma1       = gui.sigma1;        % vaccine failure (0 = perfect protection)
sigma2       = gui.sigma2;        % vaccine failure (0 = perfect protection)
vax_imm_year = gui.vax_imm_year;  % immunity duration in years

psi1 = 1 - (1 - psi1_ann)^(1/52); % weekly infant vaccination probability
psi2 = 1 - (1 - psi2_ann)^(1/52); % weekly non-infant booster probability

omega_rate = 1/(vax_imm_year*52); % waning rate of vaccine-induced immunity (weeks^-1)
omega      = 1 - exp(-omega_rate);% weekly probability of losing vaccine immunity

%% --- Initial conditions ---
% State vector format:
%   x0 = [S1; V1; I1; S2; V2; I2]
% Recovered are implicit:
%   R1 = N1 - (S1+V1+I1),    R2 = N2 - (S2+V2+I2)
x0 = gui.x0(:);



%% --- Simulation of 5 years (weekly steps), with optional inter-annual noise ---

% Pack parameters into a struct 
par = struct('N1',N1,'N2',N2,'Lambda',Lambda,'r1',r1,'r2',r2,'r3',r3,'eta',eta, ...
             'p',p,'psi1',psi1,'psi2',psi2,'sigma1',sigma1,'sigma2',sigma2, ...
             'omega',omega,'gamma1',gamma1,'gamma2',gamma2,'nu',nu, ...
             'beta11',beta11,'beta12',beta12,'beta21',beta21,'beta22',beta22, ...
             'd', d);

% Compute R0, Rc and the DFE (requires function RC_calc)
[S1_DFE, V1_DFE, S2_DFE, V2_DFE, R0, RC] = RC_calc(par); 


years       = 5;
T           = years*52;     % total number of weekly steps
year_vector = 1:years;

% Preallocation (full-horizon storage)
X          = zeros(6, T+1); % solution trajectory across all years (states include endpoints)
G1_all     = zeros(1,T);    % infection/contagion probability in age class 1 (weekly)
G2_all     = zeros(1,T);    % infection/contagion probability in age class 2 (weekly)
inc1_all   = zeros(1,T);    % incidence per week in age class 1
inc2_all   = zeros(1,T);    % incidence per week in age class 2
vax1_all   = zeros(1,T);    % vaccinations per week in age class 1
vax2_all   = zeros(1,T);    % boosters per week in age class 2
deaths_all = zeros(1,T);    % pertussis-induced deaths per week in age class 1

% Read GUI input and set noise amplitude (0 if disabled, 5% if enabled)
noise_level = 0;
if gui.noiseOn
    noise_level = 0.05;
end

% Run the model year-by-year (to allow a state perturbation at the boundary of each year)
t_counter = 0;
for yr = 1:years
    t0 = t_counter;
    t1 = t_counter + 52;

    % Integrate one year (weekly steps) starting from x0
    [x, G1, G2, inc1, inc2, vax1, vax2, deaths, R1, R2] = discrete_system(t0, t1, x0, par);

    % Place the yearly outputs into the full-horizon arrays
    X(:, t0+1:t1+1)       = x;   
    G1_all(1, t0+1:t1)    = G1; 
    G2_all(1, t0+1:t1)    = G2;
    inc1_all(1, t0+1:t1)  = inc1;
    inc2_all(1, t0+1:t1)  = inc2;
    vax1_all(1, t0+1:t1)  = vax1;
    vax2_all(1, t0+1:t1)  = vax2;

    % Prepare initial condition for next year: end-of-year state ± noise
    if yr < years
        x_end = x(:, end);
        x0    = add_noise_state(x_end, N1, N2, noise_level);
        X(:, t1+1) = x0;   
    end

    t_counter = t1;
end

%% --- Annual relevant quantities (sum weekly quantities within each year) ---
inc1_annual = zeros(years,1);
inc2_annual = zeros(years,1);
vax1_annual = zeros(years,1);
vax2_annual = zeros(years,1);

for a = 1:years
    k1 = (a-1)*52 + 1;  
    k2 = a*52;         

    inc1_annual(a) = sum(inc1_all(k1:k2));
    inc2_annual(a) = sum(inc2_all(k1:k2));

    v1 = sum(vax1_all(k1:k2));
    v2 = sum(vax2_all(k1:k2));

    vax1_annual(a) = (v1/N1)*100;    % annual percentage (infants)
    vax2_annual(a) = (v2/N2)*10000;  % annual boosters per 10,000 (non-infants)
end

%% --- Pack outputs for the GUI ---
out = struct();
out.RC          = RC;
out.R0          = R0;
out.year_vector = year_vector;
out.inc1_annual = inc1_annual;
out.inc2_annual = inc2_annual;
out.vax1_annual = vax1_annual;
out.vax2_annual = vax2_annual;
end