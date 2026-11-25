%% PerTexP simulator with optional inter-annual state noise
clear; clc; close all;

%% Optional inter-annual state noise
use_noise = false;   % true: perturb the state before each subsequent year
noise_level = 0.05;
noise_seed = [];   

% Validate optional noise settings
validateattributes(use_noise, {'logical'}, {'scalar'});
validateattributes(noise_level, {'numeric'}, ...
    {'real','scalar','finite','>=',0,'<=',1/3});
if ~isempty(noise_seed)
    validateattributes(noise_seed, {'numeric'}, ...
        {'real','scalar','integer','finite','>=',0,'<=',2^32-1});
end
if use_noise && noise_level > 0 && exist('add_noise_state','file') ~= 2
    error('pertexp_simulator:MissingNoiseFunction', ...
        'Place add_noise_state.m on the MATLAB path.');
end

%% Parameters and initial conditions
% 2024 Italy age-specific REPORTED cases (ECDC Atlas)
target1_reported = 773;
target2_reported = 3878;

% Reporting fractions
rho1 = 1;      % infants:      no under-reporting assumed
rho2 = 1/100;   % non-infants:  100x under-reporting 

% Demographic / vaccination parameters
par_base = build_par_base();

% Initial conditions
[x0, par_base, target1_true, target2_true] = build_initial_condition_and_par( ...
    par_base, rho1, rho2, target1_reported, target2_reported);

% Weekly contact rates from Mistry et al.
contact_rates = compute_contact_rates();
C11 = contact_rates.C11_week;
C12 = contact_rates.C12_week;
C21_mistry = contact_rates.C21_week;
C21 = (par_base.N1/par_base.N2)*C12; % 2024 contact reciprocity
C22 = contact_rates.C22_week;

%% Two-parameter calibration: rho_trans and q_trans

fprintf('\n--- Two-parameter transmission calibration: rho_trans + q_trans ---\n');

rho_trans_bounds = [0.01, 0.999];
q_trans_bounds   = [1e-8, 1];
weights_fit      = [1, 1];

[rho_trans, q_trans, beta11, beta12, beta21, beta22, ...
 CI1_fit, CI2_fit, err1_rel, err2_rel, ...
 slack1, slack2, J_fit, ...
 R0_final, RC_final, calib_info] = ...
    calibrate_rho_q( ...
        par_base, contact_rates, ...
        x0, target1_true, target2_true, ...
        rho_trans_bounds, q_trans_bounds, weights_fit);

fprintf('\nCalibration mode: %s\n', calib_info.mode);
fprintf('Fallback method: %s\n', calib_info.fallback_method);

fprintf('\nEstimated transmission parameters:\n');
fprintf('  rho_trans = %.8f\n', rho_trans);
fprintf('  q_trans   = %.8f per contact\n', q_trans);
fprintf('  beta11    = %.8f week^-1\n', beta11);
fprintf('  beta12    = %.8f week^-1\n', beta12);
fprintf('  beta21    = %.8f week^-1\n', beta21);
fprintf('  beta22    = %.8f week^-1\n', beta22);
fprintf('  J         = %.8e\n', J_fit);

fprintf('\nReproduction numbers:\n');
fprintf('  R0 = %.4f\n', R0_final);
fprintf('  Rc = %.4f\n', RC_final);

% Hard safety check for fallback mode.
if ~calib_info.exact_match_found
    if CI1_fit < target1_true || CI2_fit < target2_true
        error('Fallback calibration violates a required incidence lower bound.');
    end
end

par = par_base;
par.beta11 = beta11;
par.beta12 = beta12;
par.beta21 = beta21;
par.beta22 = beta22;

%% Endemic equilibrium
if RC_final > 1
    [x_endemic, endemic_info] = compute_endemic_equilibrium(par, x0);

    fprintf('\n--- Endemic equilibrium (weekly discrete system) ---\n');
    fprintf('Infants:      S1* = %.3f, V1* = %.3f, I1* = %.3f, R1* = %.3f\n', ...
        x_endemic(1), x_endemic(2), x_endemic(3), endemic_info.R1);
    fprintf('Non-infants:  S2* = %.3f, V2* = %.3f, I2* = %.3f, R2* = %.3f\n', ...
        x_endemic(4), x_endemic(5), x_endemic(6), endemic_info.R2);
    fprintf('Annualized TRUE incidence at equilibrium: infants = %.3f, non-infants = %.3f\n', ...
        endemic_info.annual_incidence(1), endemic_info.annual_incidence(2));
    fprintf('Annualized REPORTED incidence at equilibrium: infants = %.3f, non-infants = %.3f\n', ...
        rho1*endemic_info.annual_incidence(1), rho2*endemic_info.annual_incidence(2));
    fprintf('Jacobian spectral radius = %.8f -> locally stable: %d\n', ...
        endemic_info.spectral_radius, endemic_info.locally_stable);
else
    x_endemic = nan(6,1);
    endemic_info = struct();

    fprintf(['\nRc <= 1: the numerical search for an endemic ', ...
             'equilibrium was not performed.\n']);
end

%% Projection with calibrated rates and optional inter-annual state noise
years = 5;
T = years*52;

noise_active = use_noise && noise_level > 0;
if noise_active && ~isempty(noise_seed)
    rng(noise_seed, 'twister');
end
noise_info = struct('enabled',use_noise,'active',noise_active, ...
    'level',noise_level,'seed',noise_seed,'rng_before',rng, ...
    'boundary_weeks',52*(1:years-1), ...
    'state_before',zeros(6,max(years-1,0)), ...
    'state_after',zeros(6,max(years-1,0)));
fprintf('\nProjection: inter-annual noise active = %d, SD = %.3f.\n', ...
    noise_active, noise_level);

X = zeros(6, T+1); X(:,1) = x0;
inc1_all = zeros(1,T); inc2_all = zeros(1,T);
vax1_all = zeros(1,T); vax2_all = zeros(1,T); deaths_all = zeros(1,T);

x_curr = x0; t_counter = 0;
for yr = 1:years
    t0 = t_counter; t1 = t_counter + 52;
    [x, ~, ~, inc1, inc2, vax1, vax2, deaths, ~, ~] = discrete_system(t0, t1, x_curr, par);
    X(:, t0+1:t1+1)        = x;
    inc1_all(1, t0+1:t1)   = inc1;
    inc2_all(1, t0+1:t1)   = inc2;
    vax1_all(1, t0+1:t1)   = vax1;
    vax2_all(1, t0+1:t1)   = vax2;
    deaths_all(1, t0+1:t1) = deaths;
    x_curr = x(:, end);
    if yr < years
        noise_info.state_before(:,yr) = x_curr;
        if noise_active
            x_curr = add_noise_state(x_curr, par.N1, par.N2, noise_level);
        end
        noise_info.state_after(:,yr) = x_curr;
        X(:,t1+1) = x_curr;
    end
    t_counter = t1;
end

noise_info.rng_after = rng;
R1_all = par.N1 - sum(X(1:3,:),1);
R2_all = par.N2 - sum(X(4:6,:),1);

inc1_annual = zeros(years,1); inc2_annual = zeros(years,1);
vax1_annual = zeros(years,1); vax2_annual = zeros(years,1);
deaths_annual = zeros(years,1);
deaths_per1000 = nan(years,1);
for a = 1:years
    k1 = (a-1)*52+1; k2 = a*52;
    inc1_annual(a) = sum(inc1_all(k1:k2));
    inc2_annual(a) = sum(inc2_all(k1:k2));

    v1 = sum(vax1_all(k1:k2));
    v2 = sum(vax2_all(k1:k2));
    vax1_annual(a) = min((v1/par_base.N1)*100, 100);
    vax2_annual(a) = min((v2/par_base.N2)*10000, 10000);

    deaths_annual(a) = sum(deaths_all(k1:k2));
    deaths_per1000(a) = NaN;
    if inc1_annual(a) > 0
        deaths_per1000(a) = 1000 * deaths_annual(a) / inc1_annual(a);
    end
end

% Total cumulative infections over the simulation horizon
TC1 = sum(inc1_annual);
TC2 = sum(inc2_annual);
TC  = TC1 + TC2;

fprintf('\n--- Cumulative incidence in %d years ---\n', years);
fprintf('Infants:     %.0f\n', TC1);
fprintf('Non-infants: %.0f\n', TC2);
fprintf('Total:      %.0f\n', TC);



%% Plot and export the simulation results
pertexp_plots(years, inc1_annual, inc2_annual, vax1_annual, vax2_annual);


%%  Local functions
%  ------------------------------------------------------------------
function par_base = build_par_base()
    N  = 58971230;
    N1 = 380630;
    N2 = N - N1;

    eta_ann = 1/52;
    eta     = 1 - exp(-eta_ann);

    mu1_ann_per_thousand = 2.57230;
    mu1_ann              = mu1_ann_per_thousand/1000;
    r1_ann               = 1 - mu1_ann;
    r1                   = r1_ann^(1/52);

    Lambda = (1 - r1*(1 - eta)) * N1;
  
    deaths_total_2024     = 653109;
    deaths_infant_2024    = 932;
    deaths_noninfant_2024 = deaths_total_2024 - deaths_infant_2024;
    
    m2_ann = deaths_noninfant_2024 / N2;
    mu2_ann = 1 - exp(-m2_ann);
    r2_ann  = 1 - mu2_ann;
    r2 = exp(-m2_ann / 52);
    
    % Net migration inflow required to keep N2 constant (persons/week)
    M = (1 - r2) * N2 - r1 * eta * N1;
    
    if M < 0
        warning('Computed net migration inflow M is negative: %.3f persons/week.', M);
    end

    p        = 0.626;
    psi1_ann = 0.947;
    psi2_ann = 0.165e-2;
    psi1     = 1 - (1 - psi1_ann)^(1/52);
    psi2     = 1 - (1 - psi2_ann)^(1/52);

    sigma1   = 0.17;
    sigma2   = 0.08;

    omega_rate = 1/(10*52);
    omega      = 1 - exp(-omega_rate);
    nu_rate    = 1/(15*52);
    nu         = 1 - exp(-nu_rate);

    gamma1_rate = 1/3;
    gamma1      = 1 - exp(-gamma1_rate);
    gamma2_rate = 1/3;
    gamma2      = 1 - exp(-gamma2_rate);

    p_CFR = 6.3/1000;
    r3 = (r1 - p_CFR) /((1 - p_CFR) + p_CFR*(gamma1 + eta))

    d = r1 - r3;

    par_base = struct( ...
        'N1', N1, 'N2', N2, 'Lambda', Lambda, ...
        'r1', r1, 'r2', r2, 'r3', r3, 'eta', eta, 'd', d, ...
        'M', M, 'mu2_ann', mu2_ann, 'r2_ann', r2_ann, ...
        'p', p, 'psi1', psi1, 'psi2', psi2, 'sigma1', sigma1, 'sigma2', sigma2, ...
        'omega', omega, 'gamma1', gamma1, 'gamma2', gamma2, 'nu', nu, ...
        'beta11', NaN, 'beta12', NaN, 'beta21', NaN, 'beta22', NaN, ...
        'psi1_ann', psi1_ann, ...
        'gamma1_rate', gamma1_rate, 'gamma2_rate', gamma2_rate);
end

function [x0, par_base, target1_true, target2_true] = build_initial_condition_and_par( ...
        par_base, rho1, rho2, target1_reported, target2_reported)

    N1 = par_base.N1; N2 = par_base.N2;
    psi1_ann = par_base.psi1_ann; p = par_base.p;

    target1_true = target1_reported / rho1;
    target2_true = target2_reported / rho2;

    I10 = (target1_true/52) * (1/par_base.gamma1_rate)
    I20 = (target2_true/52) * (1/par_base.gamma2_rate)

    V10 = N1 * ((4/12) * p + (8/12) * psi1_ann)
    R10 = 0;
    S10 = N1 - V10 - I10 - R10

    annual_births_past10y = [397193; 407572; 414581; 431500; 450440; 469364; 483663; 495818; 510387; 517826];
    annual_vaxcoverages_past10y = [94.45; 94.76; 95.14; 94.00; 94.03; 94.99; 95.07; 94.62; 93.55; 93.33] / 100;

    V20 = dot(annual_births_past10y, annual_vaxcoverages_past10y)
 
    reported_infant_2009_2023    = 1762;
    reported_noninfant_2009_2023 = 6074;

     R20 = reported_infant_2009_2023 / rho1 ...
    + reported_noninfant_2009_2023 / rho2
    S20 = N2 - V20 - I20 - R20

    if min([S10, V10, I10, S20, V20, I20, R10, R20]) < 0
        error('build_initial_condition_and_par: negative initial condition for rho1=%.4g, rho2 = %.4g (S10=%.0f, S20 = %.0f).', rho1, rho2, S10, S20);
    end

    x0 = [S10; V10; I10; S20; V20; I20];
end
