%% PerTexP simulator (without Graphical User Interface) — optional inter-annual state noise
%
% Notes:
%  - N1(t) = N1 and N2(t) = N2 are assumed constant (infants / non-infants).
%  - Recovered compartments are NOT explicit states, but are reconstructed as:
%       R1(t) = N1 - S1(t) - V1(t) - I1(t),
%       R2(t) = N2 - S2(t) - V2(t) - I2(t).
%
% This script:
%   1) Sets model parameters and initial conditions
%   2) Runs the simulation using "discrete_system" (weekly time step)
%   3) Computes annual total quantities (cases, vaccinations, deaths)
%   4) Produces figures and saves each one as BOTH .eps and .png into two separate folders.
%
% Noise option (inter-annual variability):
%  - In the "Simulation" section you can enable/disable state noise and set its magnitude.
%  - At the end of each simulated year (except the last), the final state vector is perturbed
%    by multiplicative noise via add_noise_state() and used as the initial condition for the next year.

clear; clc; close all;

%% Plotting configuration (global defaults)
% These settings affect all figures created after this point.
set(0, 'defaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 25);
set(0, 'defaultLineMarkerSize', 2);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultFigureColor', 'w');

% Pastel colors (used consistently across plots)
c1 = [0.62, 0.80, 0.94];   % light blue  (infants)
c2 = [0.99, 0.70, 0.40];   % orange      (non-infants)

% Default figure size for single-panel bar plots
W = 700; H = 800;

%% Model parameters
% --- Demography / population sizes (Italy, January 2024) ---
N  = 58971230;              % total population
N1 = 380630;                % infants (< 1 year)
N2 = N - N1;                % non-infants (>= 1 year)

% --- Ageing from infants to non-infants (weekly probability) ---
eta_ann = 1/52;             % annual rate (1 year^-1)
eta     = 1 - exp(-eta_ann);% weekly probability of aging out of infants

% --- Survival probabilities, infants ---
mu1_ann_per_thousand = 2.57230;                   % annual death probability per 1000 live births
mu1_ann              = mu1_ann_per_thousand/1000; % annual death probability
r1_ann               = 1 - mu1_ann;               % annual survival probability
r1                   = r1_ann^(1/52);             % weekly survival probability

% --- Calibration to match desired equilibria ---
Lambda = (1 - r1*(1 - eta)) * N1;                 % weekly births into infants
r2     = 1 - (r1 * eta * N1) / N2;                % weekly survival probability, non-infants
mu2_ann = 1 - r2^52;                              % annual death probability, non-infants

% --- Vaccination (all annual inputs converted to weekly probabilities) ---
p        = 0.626;                                 % maternal vaccination coverage
psi1_ann = 0.947;                                 % annual vaccination probability for susceptible infants
psi2_ann = 0.165e-2;                              % annual booster probability for susceptible non-infants

psi1     = 1 - (1 - psi1_ann)^(1/52);             % weekly infant vaccination probability
psi2     = 1 - (1 - psi2_ann)^(1/52);             % weekly non-infant booster probability

sigma1   = 0.17;                                  % vaccine failure (infants): 0=perfect, 1=no protection
sigma2   = 0.08;                                  % vaccine failure (non-infants): 0=perfect, 1=no protection

omega_rate = 1/(10*52);                           % waning rate of vaccine-induced immunity (weeks^-1)
omega      = 1 - exp(-omega_rate);                % weekly probability of losing vaccine immunity

nu_rate = 1/(10*52);                              % waning rate of natural immunity (weeks^-1)
nu      = 1 - exp(-nu_rate);                      % weekly probability of losing natural immunity

% --- Recovery (weekly probabilities) ---
gamma1_rate = 1/3;                                % recovery rate (infants)
gamma1      = 1 - exp(-gamma1_rate);              % weekly recovery probability (infants)

gamma2_rate = 1/3;                                % recovery rate (non-infants)
gamma2      = 1 - exp(-gamma2_rate);              % weekly recovery probability (non-infants)

% --- Disease-induced mortality, infants ---
p_CFR  = 6.3/1000;                                %Case Fatality Rate (CFR)
r3 = (1 - p_CFR)/((1 - p_CFR) + p_CFR*(gamma1 + eta)); % reduced survival rate for infected infants
if r3 <= 0 || r3 >= 1
    warning('Computed r3 is outside (0,1): r3 = %.6f', r3);
end
d = r1-r3;          % weekly disease-induced death probability

% --- Transmission rates --- 
beta11 = 0.81;   beta12 = 3.0;
beta21 = 0.002;  beta22 = 0.3;

% Pack parameters into a struct
par = struct( ...
    'N1', N1, 'N2', N2, 'Lambda', Lambda, 'r1', r1, 'r2', r2, 'r3', r3, 'eta', eta, 'd', d, ...
    'p', p, 'psi1', psi1, 'psi2', psi2, 'sigma1', sigma1, 'sigma2', sigma2, ...
    'omega', omega, 'gamma1', gamma1, 'gamma2', gamma2, 'nu', nu, ...
    'beta11', beta11, 'beta12', beta12, 'beta21', beta21, 'beta22', beta22);

% Compute R0, Rc and the DFE
[S1_DFE, V1_DFE, S2_DFE, V2_DFE, R0, RC] = RC_calc(par);
fprintf('RC = %.3f,   R0 = %.3f\n', RC, R0);

%% Initial conditions
% Note: Recovered are not explicit in x0; they are computed inside/outside as:
%   R1 = N1 - S1 - V1 - I1,  R2 = N2 - S2 - V2 - I2.

% --- Infants (see Section "Parametrisation") ---
I10 = 5;
V10 = N1 * ((8/12) * psi1_ann + (2/12) * p + (2/12) * p * psi1_ann);
R10 = 0;
S10 = N1 - V10 - I10 - R10;

% --- Non-infants (see Section "Parametrisation") ---
annual_births_past10y        = [397193; 407572; 414581; 431500; 450440; 469364; 483663; 495818; 510387; 517826];
annual_vaxcoverages_past10y  = [94.76; 95.14; 94; 94.03; 94.99; 95.07; 94.62; 93.55; 93.33; 94.64] / 100;

I20 = 15;
V20 = dot(annual_births_past10y, annual_vaxcoverages_past10y);
R20 = 4064;
S20 = N2 - V20 - I20 - R20;

% Sanity check (no negative initial values)
if min([S10, V10, I10, S20, V20, I20, R10, R20]) < -1e-12
    warning('Some initial components are negative.');
end

% Initial state vector: [S1; V1; I1; S2; V2; I2]
x0 = [S10; V10; I10; S20; V20; I20];

% --- Time horizon (weekly steps) ---
t0          = 0;            % initial time (weeks)
years       = 5;            % simulation horizon (years)
T           = years * 52;   % number of weekly steps
time        = 0:T;          % time grid (weeks) [not used below, but kept for reference]
year_vector = 1:years;      % year index for bar charts

%% Preallocation (full-horizon storage)
% We simulate year-by-year, but store outputs on a single 0..T timeline.
X          = zeros(6, T+1);    % state trajectory over the full horizon (all years)
G1_all     = zeros(1, T);      % infection probability in age class 1 (weekly)
G2_all     = zeros(1, T);      % infection probability in age class 2 (weekly)
inc1_all   = zeros(1, T);      % incidence in age class 1 (weekly)
inc2_all   = zeros(1, T);      % incidence in age class 2 (weekly)
vax1_all   = zeros(1, T);      % new vaccinations in age class 1 (weekly)
vax2_all   = zeros(1, T);      % new boosters in age class 2 (weekly)
deaths_all = zeros(1, T);      % pertussis-induced deaths in infants (weekly)

%% Simulation and setting noise level
% Here you can choose to add noise or not. 
use_noise = false;   % true = add noise, false = no noise
if use_noise
    noise_level = 0.05;   % ±5%
else
    noise_level = 0;
end

% NOTE: results will differ run-to-run unless you fix RNG (e.g., rng(1)).
% This script does not set the RNG seed.

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
    deaths_all(1, t0+1:t1) = deaths;

    % Prepare initial condition for next year: end-of-year state ± noise
    if yr < years
        x_end = x(:, end);
        x0    = add_noise_state(x_end, N1, N2, noise_level);
        X(:, t1+1) = x0; 
    end

    t_counter = t1;
end

%% Annual aggregates
% New annual cases (sum of weekly incidences)
inc1_annual = zeros(years,1);
inc2_annual = zeros(years,1);

for a = 1:years
    k1 = (a-1)*52 + 1;   % first week index for year a
    k2 = a*52;           % last  week index for year a
    inc1_annual(a) = sum(inc1_all(k1:k2));
    inc2_annual(a) = sum(inc2_all(k1:k2));
end

% Annual vaccinations (sum weekly vaccinations within the year)
% Then rescale: infants in %, non-infants per 10,000.
vax1_annual = zeros(years,1);
vax2_annual = zeros(years,1);

for a = 1:years
    k1 = (a-1)*52 + 1;
    k2 = a*52;

    vax1_annual(a) = sum(vax1_all(k1:k2));
    vax2_annual(a) = sum(vax2_all(k1:k2));

    vax1_annual(a) = min((vax1_annual(a)/N1)*100, 100);
    vax2_annual(a) = min((vax2_annual(a)/N2)*10000, 10000);
end

% Annual pertussis-induced deaths (expressed per 1000 infant cases)
deaths_annual = zeros(years,1);

for a = 1:years
    k1 = (a-1)*52 + 1;
    k2 = a*52;

    deaths_annual(a) = sum(deaths_all(k1:k2));
    deaths_annual(a) = (deaths_annual(a)/inc1_annual(a))*1000;
end

%% Plots (each saved as .eps and .png)
% Output folders are created if missing.
outdir_eps = fullfile(pwd, 'fig_eps');
outdir_png = fullfile(pwd, 'fig_png');

if ~exist(outdir_eps, 'dir'); mkdir(outdir_eps); end
if ~exist(outdir_png, 'dir'); mkdir(outdir_png); end

%% 1) Annual cumulative cases by age (grouped bars)
fig = figure('Name','Annual cumulative cases by age class', 'Position',[200 200 1200 600]);
ax  = axes(fig); hold(ax,'on');

Y_inc  = [inc1_annual, inc2_annual];
hb_inc = bar(ax, year_vector, Y_inc, 'grouped');

xticks(ax, year_vector);
xticklabels(ax, year_vector);     % explicit, for clarity
grid(ax,'off'); box(ax,'off');

xlabel(ax,'Year');
ylabel(ax,'Cumulative cases');
%title(ax,'Annual cumulative cases by age class');

% Styling
hb_inc(1).FaceColor = c1;   hb_inc(1).EdgeColor = 'none'; hb_inc(1).BarWidth = 0.8;
hb_inc(2).FaceColor = c2;   hb_inc(2).EdgeColor = 'none'; hb_inc(2).BarWidth = 0.8;

legend(ax, {'Infants','Non--infants'}, 'Location','northwest', 'Interpreter','latex');

% Value labels on bars (rounded to nearest integer)
for i = 1:numel(hb_inc)
    x   = hb_inc(i).XEndPoints;
    y   = hb_inc(i).YEndPoints;
    lbl = compose('$%d$', round(y));
    text(x, y, lbl, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'Interpreter','latex', 'FontSize',18);
end

exportgraphics(fig, fullfile(outdir_eps, 'cumulative_cases_byage.eps'), 'ContentType', 'vector');
exportgraphics(ax,  fullfile(outdir_png, 'cumulative_cases_byage.png'), 'Resolution', 300);

%% 2) Annual cumulative infant cases (single bar plot)
fig = figure('Name', 'Annual cumulative infant cases', 'Position', [200, 200, W, H]);
ax  = axes(fig); hold(ax, 'on');

b = bar(ax, year_vector, inc1_annual);
b.FaceColor = c1; b.EdgeColor = 'none'; b.FaceAlpha = 1; b.BarWidth = 0.7;

xlim(ax, [0.5, 5.5]);
set(ax, 'XTick', 1:years, 'XTickLabel', 1:years);

xlabel(ax, 'Year');
ylabel(ax, 'Cumulative infant cases');

% Add headroom on y-axis to avoid labels being clipped
yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

% Value labels above bars
xt  = b.XEndPoints;
yt  = b.YEndPoints;
dy  = 0.005 * (yl(2) - yl(1));
lbl = compose('$%d$', round(b.YData));

text(xt, yt + dy, lbl, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize',18, 'Interpreter','latex', 'Clipping','off');

exportgraphics(fig, fullfile(outdir_eps, 'cumulative_cases_infants.eps'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir_png, 'cumulative_cases_infants.png'), 'Resolution', 300);

%% 3) Annual cumulative non-infant cases (single bar plot)
fig = figure('Name', 'Annual cumulative non-infant cases', 'Position', [200, 200, W, H]);
ax  = axes(fig); hold(ax, 'on');

b = bar(ax, year_vector, inc2_annual);
b.FaceColor = c2; b.EdgeColor = 'none'; b.FaceAlpha = 1; b.BarWidth = 0.7;

xlim(ax, [0.5, 5.5]);
set(ax, 'XTick', 1:years, 'XTickLabel', 1:years);

xlabel(ax, 'Year');
ylabel(ax, 'Cumulative non--infant cases');

yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

xt  = b.XEndPoints;
yt  = b.YEndPoints;
dy  = 0.005 * (yl(2) - yl(1));
lbl = compose('$%d$', round(b.YData));

text(xt, yt + dy, lbl, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize',18, 'Interpreter','latex', 'Clipping','off');

exportgraphics(fig, fullfile(outdir_eps, 'cumulative_cases_noninfants.eps'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir_png, 'cumulative_cases_noninfants.png'), 'Resolution', 300);

%% 4) Infant vaccinations (annual %, single bar plot)
fig = figure('Name', 'Infant vaccinations (annual %)', 'Position', [200, 200, 600, 600]);
ax  = axes(fig); hold(ax, 'on');

b = bar(ax, year_vector, vax1_annual);
b.FaceColor = c1; b.EdgeColor = 'none'; b.FaceAlpha = 1; b.BarWidth = 0.8;

xlim(ax, [0.5, 5.5]);
set(ax, 'XTick', 1:years, 'XTickLabel', 1:years);

xlabel(ax, 'Year');
ylabel(ax, 'Infant vaccinations (\%)');

yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

xt  = b.XEndPoints;
yt  = b.YEndPoints;
dy  = 0.005 * (yl(2) - yl(1));
lbl = compose('$%.2f$', b.YData);

text(xt, yt + dy, lbl, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize',23, 'Interpreter','latex', 'Clipping','off');

exportgraphics(ax, fullfile(outdir_eps, 'vaccinations_infants_percent.eps'), 'ContentType', 'vector');
exportgraphics(ax, fullfile(outdir_png, 'vaccinations_infants_percent.png'), 'Resolution', 300);

%% 5) Non-infant boosters (per 10000, single bar plot)
fig = figure('Name', 'Non-infant boosters (per 10,000)', 'Position', [200, 200, 600, 600]);
ax  = axes(fig); hold(ax, 'on');

b = bar(ax, year_vector, vax2_annual);
b.FaceColor = c2; b.EdgeColor = 'none'; b.FaceAlpha = 1; b.BarWidth = 0.8;

xlim(ax, [0.5, 5.5]);
set(ax, 'XTick', 1:years, 'XTickLabel', 1:years);

xlabel(ax, 'Year');
ylabel(ax, 'Non-infant boosters (per 10\,000)');

yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

xt  = b.XEndPoints;
yt  = b.YEndPoints;
dy  = 0.005 * (yl(2) - yl(1));
lbl = compose('$%.2f$', b.YData);

text(xt, yt + dy, lbl, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize',23, 'Interpreter','latex', 'Clipping','off');

exportgraphics(fig, fullfile(outdir_eps, 'boosters_noninfants_per10000.eps'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir_png, 'boosters_noninfants_per10000.png'), 'Resolution', 300);



%% 2) Annual cumulative pertussis-related deaths 
fig = figure('Name', 'Annual cumulative pertussis-related deaths', 'Position', [200, 200, W, H]);
ax  = axes(fig); hold(ax, 'on');

b = bar(ax, year_vector, deaths_annual);
b.FaceColor = c1; b.EdgeColor = 'none'; b.FaceAlpha = 1; b.BarWidth = 0.7;

xlim(ax, [0.5, 5.5]);
set(ax, 'XTick', 1:years, 'XTickLabel', 1:years);

xlabel(ax, 'Year');
ylabel(ax, 'Cumulative pertussis-related deaths (per 1000)');

% Add headroom on y-axis to avoid labels being clipped
yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

% Value labels above bars
xt  = b.XEndPoints;
yt  = b.YEndPoints;
dy  = 0.005 * (yl(2) - yl(1));
lbl = compose('$%d$', round(b.YData));

text(xt, yt + dy, lbl, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize',18, 'Interpreter','latex', 'Clipping','off');

exportgraphics(fig, fullfile(outdir_eps, 'cumulative_pertussis_deaths.eps'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir_png, 'cumulative_pertussis_deaths.png'), 'Resolution', 300);

%% Summary
disp(['EPS saved in: ', outdir_eps]);
disp(['PNG saved in: ', outdir_png]);
