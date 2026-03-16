% Contour plots of Rc
%
%% 1) As a function of (psi1_ann, p), with all other parameters fixed
clear; clc; close all;

% Figure size and colors (defaults)
set(0,'defaultLineLineWidth', 2)
set(0,'defaultLineColor', 'b')
set(0,'defaultAxesFontSize',25)
set(0,'defaultLineMarkerSize',2)

% Typography (LaTeX rendering)
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot, 'defaultFigureColor', 'w');

s  = 40;                  % marker size (pt^2 for scatter)

% Pastel colors
c1 = [0.62  0.80 0.94];   % light blue
c2 = [0.99 0.70 0.40];    % orange


% Parameters

% Target population (Italy, January 2024)
N  = 58971230;    % total population
N1 = 380630;      % < 1 year
N2 = N-N1;        % >= 1 year

% Fraction of infants who turn 1 year old per week
eta_ann = 1/52;
eta     = 1-exp(-eta_ann);

% Survival probabilities
mu1_ann_per_thousand = 2.57230;
mu1_ann              = mu1_ann_per_thousand/1000; 
r1_ann               = 1 - mu1_ann; 
r1                   = r1_ann^(1/52);        

% Calibration to match the desired equilibria
Lambda = (1 - r1*(1-eta)) * N1;          
r2     = 1 - (r1*eta*N1) / N2;                  
mu2_ann = 1-r2^52;                         

% Vaccination
p        = 0.626;                              
psi1_ann = 0.947;                         
psi2_ann = 0.165e-2;                             
psi1     = 1-(1-psi1_ann)^(1/52);             
psi2     = 1-(1-psi2_ann)^(1/52);            
sigma1   = 0.17;                           
sigma2   = 0.08;                              
omega_rate = 1/(10*52);                  
omega      = 1-exp(-omega_rate);          
nu_rate    = 1/(10*52);        
nu         = 1-exp(-nu_rate);

% Recovery
gamma1_rate = 1/3;              
gamma1      = 1-exp(-gamma1_rate);       
gamma2_rate = 1/3;          
gamma2      = 1-exp(-gamma2_rate);

% --- Disease-induced mortality, infants ---
p_CFR  = 6.3/1000;                                     % Case Fatality Rate (CFR)
r3 = (1 - p_CFR)/((1 - p_CFR) + p_CFR*(gamma1 + eta)); % reduced survival rate for infected infants
if r3 <= 0 || r3 >= 1
    warning('Computed r3 is outside (0,1): r3 = %.6f', r3);
end
d = r1-r3;          % weekly disease-induced death probability

% Transmission rates
beta11 = 0.81;    beta12 = 3.0;      
beta21 = 0.002;   beta22 = 0.3;    

% Pack parameters into a struct
par = struct('N1',N1,'N2',N2,'Lambda',Lambda,'r1',r1,'r2',r2,'r3',r3,'eta',eta, ...
             'p',p,'psi1',psi1,'psi2',psi2,'sigma1',sigma1,'sigma2',sigma2, ...
             'omega',omega,'gamma1',gamma1,'gamma2',gamma2,'nu',nu, ...
             'beta11',beta11,'beta12',beta12,'beta21',beta21,'beta22',beta22, ...
             'd', d);

% Save plots
outdir = fullfile(pwd,'fig_contour_eps');
if ~exist(outdir,'dir'), mkdir(outdir); end


% --- Grids for (psi1_ann, p) ---
psi1_ann_grid = linspace(0,1);   % annual infant coverage (baseline ~0.947) [
p_grid        = linspace(0,1);   % maternal coverage

% Allocate Rc matrix
RCmat = zeros(length(p_grid), length(psi1_ann_grid));

% --- Loop over the grid ---
for i = 1:length(p_grid)
    p = p_grid(i);
    for j = 1:length(psi1_ann_grid)
        psi1_ann = psi1_ann_grid(j);
        psi1 = 1 - (1 - psi1_ann)^(1/52);  % annual -> weekly

        % Pack parameters for RC_calc
        par = struct();
        par.N1 = N1; par.N2 = N2;
        par.Lambda = Lambda; par.r1 = r1; par.r2 = r2; par.r3 = r3; par.eta = eta;
        par.p = p; par.psi1 = psi1; par.psi2 = psi2;
        par.sigma1 = sigma1; par.sigma2 = sigma2; par.omega = omega;
        par.gamma1 = gamma1; par.gamma2 = gamma2; par.nu = nu;
        par.beta11 = beta11; par.beta12 = beta12; par.beta21 = beta21; par.beta22 = beta22;

        % Compute Rc
        [~, ~, ~, ~, ~, RC] = RC_calc(par);
        RCmat(i,j) = RC;
    end
end

% --- Contour plot ---
W = 700; H = 700;
fig1 = figure('Name','RC contour','Position',[200 200 W H]);
t = tiledlayout(fig1, 1, 1, 'TileSpacing','compact','Padding','compact');
ax = nexttile(t); hold(ax,'on');

[X,Y] = meshgrid(psi1_ann_grid, p_grid);   % X=psi1_ann, Y=p

% Filled contour map + isolines
contourf(ax, X, Y, RCmat, 20, 'LineStyle','none');
cb = colorbar(ax); cb.Label.Interpreter = 'latex'; cb.TickLabelInterpreter = 'latex';

% Highlight the iso-curve Rc = 1
[C1,h1] = contour(ax, X, Y, RCmat, [1 1], 'k', 'LineWidth', 3);
clabel(C1,h1,'Color','k','Interpreter','latex');

% Labels and ticks (psi1_ann shown in percent)
xlabel(ax, ' Annual infant coverage, $\psi_{1}^{\mathrm{ann}}$','Interpreter','latex','FontSize',30);
ylabel(ax, ' Annual maternal coverage, $p$','Interpreter','latex','FontSize',30);
title(' Control reproduction number, $\mathcal{R}_c$','Interpreter','latex','FontSize',22);
xt = get(ax,'XTick');
xticklabels(ax, compose('%.0f\\%%', xt*100)); % 0–100%
yt = get(ax,'YTick');
yticklabels(ax, compose('%.0f\\%%', yt*100)); % 0–100%
grid(ax,'on'); ax.Box='off'; ax.TickDir='out';

% Baseline lines
psi1_ann_base = 0.947;
p_base        = 0.626;
xline(ax, psi1_ann_base, ':', 'Color',[1 1 1], 'LineWidth',3);
yline(ax, p_base,        ':', 'Color',[1 1 1], 'LineWidth',3);

exportgraphics(ax, fullfile(outdir,'RC_contour_psi1.eps'), 'ContentType','vector');


%% 2) As a function of (psi2_ann, p), with all other parameters fixed

% Vaccination
p        = 0.626;                                % maternal vaccination coverage
psi1_ann = 0.947;                                % annual vaccination probability of susceptible infants
psi2_ann = 0.165e-2;                             % annual booster probability of susceptible non-infants
psi1     = 1-(1-psi1_ann)^(1/52);                % weekly vaccination probability for susceptible infants

% Pack parameters into a struct
par = struct('N1',N1,'N2',N2,'Lambda',Lambda,'r1',r1,'r2',r2,'r3',r3,'eta',eta, ...
             'p',p,'psi1',psi1,'psi2',psi2,'sigma1',sigma1,'sigma2',sigma2, ...
             'omega',omega,'gamma1',gamma1,'gamma2',gamma2,'nu',nu, ...
             'beta11',beta11,'beta12',beta12,'beta21',beta21,'beta22',beta22, ...
             'd', d);

% --- Grids for (psi2_ann, p) ---
psi2_ann_grid = linspace(0,0.1);   % annual non-infant booster coverage
p_grid        = linspace(0,1);   % maternal coverage

% Allocate Rc matrix
RCmat = zeros(length(p_grid), length(psi2_ann_grid));

% --- Loop over the grid ---
for i = 1:length(p_grid)
    p = p_grid(i);
    for j = 1:length(psi2_ann_grid)
        psi2_ann = psi2_ann_grid(j);
        psi2 = 1 - (1 - psi2_ann)^(1/52);

        % Pack parameters for RC_calc
        par = struct();
        par.N1 = N1; par.N2 = N2;
        par.Lambda = Lambda; par.r1 = r1; par.r2 = r2; par.r3 = r3; par.eta = eta;
        par.p = p; par.psi1 = psi1; par.psi2 = psi2;
        par.sigma1 = sigma1; par.sigma2 = sigma2; par.omega = omega;
        par.gamma1 = gamma1; par.gamma2 = gamma2; par.nu = nu;
        par.beta11 = beta11; par.beta12 = beta12; par.beta21 = beta21; par.beta22 = beta22;

        % Compute Rc
        [~, ~, ~, ~, ~, RC] = RC_calc(par);
        RCmat(i,j) = RC;
    end
end

% --- Contour plot ---
fig2 = figure('Name','RC contour booster' ,'Position',[200 200 W H]);
t = tiledlayout(fig2, 1, 1, 'TileSpacing','compact','Padding','compact');
ax = nexttile(t); hold(ax,'on');

[X,Y] = meshgrid(psi2_ann_grid, p_grid);   % X=psi2_ann, Y=p

% Filled contour map + isolines
contourf(ax, X, Y, RCmat, 20, 'LineStyle','none');
cb = colorbar(ax); cb.Label.Interpreter = 'latex';

% Highlight the iso-curve Rc = 1
[C1,h1] = contour(ax, X, Y, RCmat, [1 1], 'k', 'LineWidth', 3);
clabel(C1,h1,'Color','k','Interpreter','latex');

% Labels and ticks (percent format)
xlabel(ax, ' Annual booster coverage, $\psi_{2}^{\mathrm{ann}}$','Interpreter','latex','FontSize',30);
ylabel(ax, ' Annual maternal coverage, $p$','Interpreter','latex','FontSize',30);
title(' Control reproduction number, $\mathcal{R}_c$','Interpreter','latex','FontSize',22);
xt = get(ax,'XTick');
xticklabels(ax, compose('%.0f\\%%', xt*100)); % 0–100%
yt = get(ax,'YTick');
yticklabels(ax, compose('%.0f\\%%', yt*100)); % 0–100%
grid(ax,'on'); ax.Box='off'; ax.TickDir='out';

% Baseline lines
psi2_ann_base = 16.5/10000;   % = 0.00165
p_base        = 0.626;
xline(ax, psi2_ann_base, ':', 'Color',[1 1 1], 'LineWidth',3);
yline(ax, p_base,':', 'Color',[1 1 1], 'LineWidth',3);

exportgraphics(ax, fullfile(outdir,'RC_contour_psi2.eps'), 'ContentType','vector');


disp(['EPS saved in: ' outdir]);
