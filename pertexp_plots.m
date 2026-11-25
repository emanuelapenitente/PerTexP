function pertexp_plots(years, inc1_annual, inc2_annual, ...
    vax1_annual, vax2_annual)
%PERTEXP_PLOTS Plot and export annual cumulative infections and vaccinations.
%
% Inputs:
%   years        : simulation horizon in years.
%   inc1_annual  : annual cumulative infant infections.
%   inc2_annual  : annual cumulative non-infant infections.
%   vax1_annual  : annual infant vaccinations, expressed as percentages.
%   vax2_annual  : annual non-infant boosters, per 10,000 individuals.
%
% Generates four separate figures and exports them as EPS and PNG
% into fig_eps and fig_png in the working directory.

year_vector = 1:years;

%% Plotting configuration
set(0, 'defaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 25);
set(0, 'defaultLineMarkerSize', 2);

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesToolbarVisible', 'off');

c1 = [0.62, 0.80, 0.94];   % Infants
c2 = [0.99, 0.70, 0.40];   % Non-infants

%% Output directories
outdir_png = fullfile(pwd, 'fig_png');
outdir_eps = fullfile(pwd, 'fig_eps');

if ~exist(outdir_png, 'dir')
    mkdir(outdir_png);
end
if ~exist(outdir_eps, 'dir')
    mkdir(outdir_eps);
end

%% 1) Annual cumulative infections
incidence_data = {inc1_annual(:).', inc2_annual(:).'};
age_names = {'Infants', 'Non-infants'};
age_colours = {c1, c2};
age_files = {'cum_inc_infants', 'cum_inc_noninfants'};

for age = 1:2
    validateattributes(incidence_data{age}, {'numeric'}, ...
        {'real', 'finite', 'nonnegative', 'numel', years});
end

for age = 1:2
    fig = figure( ...
        'Name', ['Estimated cumulative infections: ' age_names{age}], ...
        'Position', [200, 200, 600, 600]);

    ax = axes(fig);
    hold(ax, 'on');

    b = bar(ax, year_vector, incidence_data{age});
    b.FaceColor = age_colours{age};
    b.EdgeColor = 'none';
    b.FaceAlpha = 1;
    b.BarWidth = 0.8;

    xlim(ax, [0.5, years + 0.5]);
    set(ax, 'XTick', year_vector, 'XTickLabel', year_vector);

    xlabel(ax, 'Year');
    ylabel(ax, 'Estimated cumulative infections');

    yl = ylim(ax);
    yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
    ylim(ax, yl);

    xt = b.XEndPoints;
    yt = b.YEndPoints;
    dy = 0.005 * (yl(2) - yl(1));
    lbl = compose('$%d$', round(b.YData));

    text(ax, xt, yt + dy, lbl, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 23, ...
        'Interpreter', 'latex', ...
        'Clipping', 'off');

    exportgraphics(ax, ...
        fullfile(outdir_eps, [age_files{age} '.eps']), ...
        'ContentType', 'vector');

    exportgraphics(ax, ...
        fullfile(outdir_png, [age_files{age} '.png']), ...
        'Resolution', 300);
end

%% 2) Annual cumulative infant vaccinations (%)
fig = figure( ...
    'Name', 'Infant vaccinations (annual %)', ...
    'Position', [200, 200, 600, 600]);

ax = axes(fig);
hold(ax, 'on');

b = bar(ax, year_vector, vax1_annual);
b.FaceColor = c1;
b.EdgeColor = 'none';
b.FaceAlpha = 1;
b.BarWidth = 0.8;

xlim(ax, [0.5, years + 0.5]);
set(ax, 'XTick', year_vector, 'XTickLabel', year_vector);

xlabel(ax, 'Year');
ylabel(ax, 'Infant vaccinations (\%)');

yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

xt = b.XEndPoints;
yt = b.YEndPoints;
dy = 0.005 * (yl(2) - yl(1));
lbl = compose('$%.2f$', b.YData);

text(ax, xt, yt + dy, lbl, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 23, ...
    'Interpreter', 'latex', ...
    'Clipping', 'off');

exportgraphics(ax, ...
    fullfile(outdir_eps, 'vax_infants_percent.eps'), ...
    'ContentType', 'vector');

exportgraphics(ax, ...
    fullfile(outdir_png, 'vax_infants_percent.png'), ...
    'Resolution', 300);

%% 3) Annual cumulative non-infant boosters (per 10,000)
fig = figure( ...
    'Name', 'Non-infant boosters (per 10,000)', ...
    'Position', [200, 200, 600, 600]);

ax = axes(fig);
hold(ax, 'on');

b = bar(ax, year_vector, vax2_annual);
b.FaceColor = c2;
b.EdgeColor = 'none';
b.FaceAlpha = 1;
b.BarWidth = 0.8;

xlim(ax, [0.5, years + 0.5]);
set(ax, 'XTick', year_vector, 'XTickLabel', year_vector);

xlabel(ax, 'Year');
ylabel(ax, 'Non-infant boosters (per 10\,000)');

yl = ylim(ax);
yl = [yl(1), max([yl(2), 1.15 * max(b.YData)])];
ylim(ax, yl);

xt = b.XEndPoints;
yt = b.YEndPoints;
dy = 0.005 * (yl(2) - yl(1));
lbl = compose('$%.2f$', b.YData);

text(ax, xt, yt + dy, lbl, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 23, ...
    'Interpreter', 'latex', ...
    'Clipping', 'off');

exportgraphics(fig, ...
    fullfile(outdir_eps, 'boosters_noninfants_per10000.eps'), ...
    'ContentType', 'vector');

exportgraphics(fig, ...
    fullfile(outdir_png, 'boosters_noninfants_per10000.png'), ...
    'Resolution', 300);

end