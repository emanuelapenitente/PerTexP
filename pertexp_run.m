function out = pertexp_run(gui)
%PERTEXP_RUN Run the PerTexP model from the graphical interface.
%   settings = pertexp_run('defaults') returns baseline GUI inputs.
%   out = pertexp_run(settings) simulates five years.
%
%   Transmission rates are FIXED at the full-precision calibrated baseline
%   values represented (rounded) in Table 2 of the main manuscript.
%   No calibration is performed by this function, including on the first run.
%   Changing vaccination parameters, initial conditions or noise never
%   changes beta11, beta12, beta21 or beta22.
%
%   gui fields: p, psi1_ann, psi2_ann, sigma1, sigma2 (fractions),
%   vax_imm_year (years), x0 = [S1;V1;I1;S2;V2;I2], noiseOn (logical).
%   Optional: noise_level (default 0.05), noise_seed (default []).
%   Recovered individuals are obtained by population balance constraint.
%   This GUI runner does not require calibration tools or contact data.

par_base = build_par_base();
[baseline_x0, par_base, ~, ~] = ...
    build_initial_condition_and_par(par_base, 1, 1/100, 773, 3878);
defaults = struct('p', par_base.p, 'psi1_ann', par_base.psi1_ann, ...
    'psi2_ann', 0.00165, 'sigma1', par_base.sigma1, ...
    'sigma2', par_base.sigma2, 'vax_imm_year', 10, ...
    'x0', baseline_x0, 'noiseOn', false, 'noise_level', 0.05, ...
    'noise_seed', [], 'N1', par_base.N1, 'N2', par_base.N2);
if nargin == 1 && (ischar(gui) || (isstring(gui) && isscalar(gui)))
    if strcmpi(gui, 'defaults')
        out = defaults;
        return
    end
    error('pertexp_run:Command', 'The supported command is ''defaults''.');
end
if nargin == 0
    gui = defaults;
end
validateattributes(gui, {'struct'}, {'scalar'});
fraction_fields = {'p','psi1_ann','psi2_ann','sigma1','sigma2'};
for k = 1:numel(fraction_fields)
    name = fraction_fields{k};
    validateattributes(gui.(name), {'numeric'}, ...
        {'real','finite','scalar','>=',0,'<=',1}, mfilename, name);
end
validateattributes(gui.vax_imm_year, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'vax_imm_year');
validateattributes(gui.x0, {'numeric'}, ...
    {'real','finite','vector','numel',6,'nonnegative'}, mfilename, 'x0');
validateattributes(gui.noiseOn, {'logical'}, {'scalar'});
x0 = gui.x0(:);
if sum(x0(1:3)) > par_base.N1 + 1e-8 || ...
        sum(x0(4:6)) > par_base.N2 + 1e-8
    error('pertexp_run:InitialPopulation', ...
        'S + V + I cannot exceed the population of the corresponding age group.');
end
if ~isfield(gui, 'noise_level'); gui.noise_level = 0.05; end
if ~isfield(gui, 'noise_seed'); gui.noise_seed = []; end
validateattributes(gui.noise_level, {'numeric'}, ...
    {'real','finite','scalar','>=',0,'<=',1/3});
if ~isempty(gui.noise_seed)
    validateattributes(gui.noise_seed, {'numeric'}, ...
        {'real','finite','scalar','integer','>=',0,'<=',2^32-1});
end
noise_active = gui.noiseOn && gui.noise_level > 0;
if noise_active && exist('add_noise_state', 'file') ~= 2
    error('pertexp_run:MissingNoiseFunction', ...
        'Place the current add_noise_state.m on the MATLAB path.');
end

%% Fixed transmission rates: full-precision baseline (units: week^-1)
% Table 2 displays these as 0.006, 0.34, 0.002 and 0.31, respectively.
% Keep the full precision here to retain the current baseline projections.
% These constants do not depend on any input from the GUI.
par = par_base;
par.beta11 = 0.006313089771692;
par.beta12 = 0.336117274809097;
par.beta21 = 0.002183563887562;
par.beta22 = 0.313677852807034;

%% Apply vaccination scenario inputs
par.p = gui.p;
par.psi1_ann = gui.psi1_ann;
par.psi2_ann = gui.psi2_ann;
par.psi1 = 1 - (1 - gui.psi1_ann)^(1/52);
par.psi2 = 1 - (1 - gui.psi2_ann)^(1/52);
par.sigma1 = gui.sigma1;
par.sigma2 = gui.sigma2;
par.omega = 1 - exp(-1/(gui.vax_imm_year*52));
[S1_DFE, V1_DFE, S2_DFE, V2_DFE, R0, RC] = RC_calc(par);

%% Five-year projection; noise acts only at the four year boundaries
years = 5;
T = years*52;
if noise_active && ~isempty(gui.noise_seed)
    rng(gui.noise_seed, 'twister');
end
noise_info = struct('enabled', gui.noiseOn, 'active', noise_active, ...
    'level', gui.noise_level, 'seed', gui.noise_seed, 'rng_before', rng, ...
    'boundary_weeks', 52*(1:years-1), ...
    'state_before', zeros(6,years-1), 'state_after', zeros(6,years-1));
X = zeros(6,T+1);
inc1_all = zeros(1,T); inc2_all = zeros(1,T);
vax1_all = zeros(1,T); vax2_all = zeros(1,T);
x_curr = x0;
for yr = 1:years
    t0 = (yr-1)*52; t1 = yr*52;
    [x, ~, ~, inc1, inc2, vax1, vax2, ~, R1, R2] = ...
        discrete_system(t0, t1, x_curr, par);
    if any(~isfinite(x(:))) || any(x(:) < -1e-7) || ...
            any(R1 < -1e-7) || any(R2 < -1e-7)
        error('pertexp_run:InvalidTrajectory', ...
            ['These inputs produce nonphysical compartment values. ', ...
             'Please revise the vaccination parameters or initial conditions.']);
    end
    X(:,t0+1:t1+1) = x;
    inc1_all(t0+1:t1) = inc1; inc2_all(t0+1:t1) = inc2;
    vax1_all(t0+1:t1) = vax1; vax2_all(t0+1:t1) = vax2;
    x_curr = x(:,end);
    if yr < years
        noise_info.state_before(:,yr) = x_curr;
        if noise_active
            x_curr = add_noise_state(x_curr, par.N1, par.N2, gui.noise_level);
        end
        noise_info.state_after(:,yr) = x_curr;
        X(:,t1+1) = x_curr;
    end
end
noise_info.rng_after = rng;

out = struct();
out.RC = RC; out.R0 = R0; out.year_vector = 1:years;
out.inc1_annual = sum(reshape(inc1_all,52,years),1).';
out.inc2_annual = sum(reshape(inc2_all,52,years),1).';
% Retain the same reporting conventions and caps as pertexp_simulator.
out.vax1_annual = min(sum(reshape(vax1_all,52,years),1).'/par.N1*100,100);
out.vax2_annual = min(sum(reshape(vax2_all,52,years),1).'/par.N2*10000,10000);
out.X = X;
out.par = par;
out.x0 = x0;
out.DFE = [S1_DFE; V1_DFE; S2_DFE; V2_DFE];
out.beta = [par.beta11 par.beta12; par.beta21 par.beta22];
out.noise_info = noise_info;
end

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
    r3 = (r1 - p_CFR) /((1 - p_CFR) + p_CFR*(gamma1 + eta));

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

    I10 = (target1_true/52) * (1/par_base.gamma1_rate);
    I20 = (target2_true/52) * (1/par_base.gamma2_rate);

    V10 = N1 * ((4/12) * p + (8/12) * psi1_ann);
    R10 = 0;
    S10 = N1 - V10 - I10 - R10;

    annual_births_past10y = [397193; 407572; 414581; 431500; 450440; 469364; 483663; 495818; 510387; 517826];
    annual_vaxcoverages_past10y = [94.45; 94.76; 95.14; 94.00; 94.03; 94.99; 95.07; 94.62; 93.55; 93.33] / 100;

    V20 = dot(annual_births_past10y, annual_vaxcoverages_past10y);
 
    reported_infant_2009_2023    = 1762;
    reported_noninfant_2009_2023 = 6074;

     R20 = reported_infant_2009_2023 / rho1 ...
    + reported_noninfant_2009_2023 / rho2;
    S20 = N2 - V20 - I20 - R20;

    if min([S10, V10, I10, S20, V20, I20, R10, R20]) < 0
        error('build_initial_condition_and_par: negative initial condition for rho1=%.4g, rho2 = %.4g (S10=%.0f, S20 = %.0f).', rho1, rho2, S10, S20);
    end

    x0 = [S10; V10; I10; S20; V20; I20];
end
