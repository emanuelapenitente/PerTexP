function results = RC_contour(par_baseline)
%RC_CONTOUR Vaccination contours using the calibrated Italian baseline.
%
% First run pertexp_simulator.m, then call:
%   results_RC = RC_contour(par);
%
% Pass the calibrated structure par, NOT par_base (whose betas are NaN).
% No calibration is performed here: the four baseline betas stay fixed.
% All other baseline parameters, including demography and waning, are copied.
%
% Figures:
%   1) Rc versus infant annual coverage psi1_ann and maternal coverage p.
%   2) Rc versus booster annual coverage psi2_ann (0--5%) and maternal p.
% The other vaccination coverage is fixed at its baseline value in each plot.
%
% Outputs are saved in fig_contour_eps and fig_contour_png.
% Both figures use the same canvas and axes positions for side-by-side TeX.

    if nargin ~= 1 || ~isstruct(par_baseline) || ~isscalar(par_baseline)
        error('RC_contour:Input', ...
            'Run the main simulator first, then call RC_contour(par).');
    end

    required_fields = {'N1','N2','Lambda','M','r1','r2','r3','eta', ...
        'p','psi1','psi2','sigma1','sigma2','omega','gamma1','gamma2', ...
        'beta11','beta12','beta21','beta22'};
    for k = 1:numel(required_fields)
        field = required_fields{k};
        if ~isfield(par_baseline, field)
            error('RC_contour:MissingField', ...
                'The baseline structure is missing field %s.', field);
        end
        validateattributes(par_baseline.(field), {'numeric'}, ...
            {'real','scalar','finite'}, mfilename, ['par.' field]);
    end
    for field = {'p','psi1','psi2'}
        validateattributes(par_baseline.(field{1}), {'numeric'}, ...
            {'>=',0,'<=',1}, mfilename, ['par.' field{1}]);
    end
    validateattributes([par_baseline.beta11, par_baseline.beta12, ...
        par_baseline.beta21, par_baseline.beta22], ...
        {'numeric'}, {'nonnegative','finite'});

    % Derive annual coverages from the probabilities actually used in par.
    p_base = par_baseline.p;
    psi1_ann_base = 1 - (1 - par_baseline.psi1)^52;
    psi2_ann_base = 1 - (1 - par_baseline.psi2)^52;

    [~,~,~,~,R0_base,RC_base] = RC_calc(par_baseline);
    fprintf('\n--- Rc contours: calibrated transmission rates held fixed ---\n');
    fprintf('beta11 = %.10g; beta12 = %.10g; beta21 = %.10g; beta22 = %.10g\n', ...
        par_baseline.beta11, par_baseline.beta12, ...
        par_baseline.beta21, par_baseline.beta22);
    fprintf('Baseline: p = %.6g; psi1_ann = %.6g; psi2_ann = %.6g\n', ...
        p_base, psi1_ann_base, psi2_ann_base);
    fprintf('Baseline reproduction numbers: R0 = %.8f; Rc = %.8f\n', ...
        R0_base, RC_base);

    % Include the baseline coordinates explicitly in both grids.
    p_grid = unique([linspace(0,1,100), p_base]);
    psi_grids = {unique([linspace(0,1,100), psi1_ann_base]), ...
                 unique([linspace(0,0.05,100), psi2_ann_base])};
    weekly_fields = {'psi1','psi2'};
    annual_fields = {'psi1_ann','psi2_ann'};
    annual_baselines = [psi1_ann_base, psi2_ann_base];
    xlabels = {{'Annual infant vaccination'; 'coverage, $\psi_1^{\mathrm{ann}}$'}, ...
               {'Annual booster vaccination'; 'coverage, $\psi_2^{\mathrm{ann}}$'}};
    names = {'RC_contour_psi1','RC_contour_psi2'};

    outdir_eps = fullfile(pwd, 'fig_contour_eps');
    outdir_png = fullfile(pwd, 'fig_contour_png');
    if ~exist(outdir_eps,'dir'), mkdir(outdir_eps); end
    if ~exist(outdir_png,'dir'), mkdir(outdir_png); end

    results = struct();
    results.par_baseline = par_baseline;
    results.R0_baseline = R0_base;
    results.RC_baseline = RC_base;
    results.p_grid = p_grid;
    results.figures = gobjects(1,2);

    for panel = 1:2
        psi_ann_grid = psi_grids{panel};
        RCmat = zeros(numel(p_grid), numel(psi_ann_grid));

        for i = 1:numel(p_grid)
            for j = 1:numel(psi_ann_grid)
                par_grid = par_baseline;
                par_grid.p = p_grid(i);
                par_grid.(weekly_fields{panel}) = ...
                    1 - (1 - psi_ann_grid(j))^(1/52);
                par_grid.(annual_fields{panel}) = psi_ann_grid(j);
                [~,~,~,~,~,RCmat(i,j)] = RC_calc(par_grid);
            end
        end

        % Check that the baseline point reproduces the simulator's Rc.
        [~,ip] = min(abs(p_grid - p_base));
        [~,jp] = min(abs(psi_ann_grid - annual_baselines(panel)));
        baseline_error = abs(RCmat(ip,jp) - RC_base);
        if baseline_error > 1e-8 * max(1,abs(RC_base))
            error('RC_contour:BaselineMismatch', ...
                'The baseline contour value does not reproduce Rc.');
        end

        fig = figure('Name', names{panel}, 'Units','pixels', ...
            'Position', [200 200 820 900], 'Color','w', ...
            'ToolBar','none','MenuBar','none','Renderer','painters');
        % Reserve identical margins for two-line labels and the colour bar.
        axes_position = [0.23 0.23 0.57 0.70];
        ax = axes('Parent',fig,'Units','normalized', ...
            'Position',axes_position,'PositionConstraint','innerposition');
        hold(ax,'on');
        [Xgrid,Ygrid] = meshgrid(psi_ann_grid,p_grid);
        contourf(ax,Xgrid,Ygrid,RCmat,20,'LineStyle','none');
        colormap(ax,parula(256));
        cb = colorbar(ax);
        cb.Units = 'normalized';
        cb.Position = [0.835 0.23 0.025 0.70];
        ax.Position = axes_position;
        cb.FontSize = 25;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.Interpreter = 'latex';
        % Put Rc above the bar to avoid a long right-hand margin.
        title(cb,'$\mathcal{R}_c$','Interpreter','latex','FontSize',25);

        if min(RCmat(:)) < 1 && max(RCmat(:)) > 1
            [C,h] = contour(ax,Xgrid,Ygrid,RCmat,[1 1], ...
                'Color','k','LineWidth',2);
            clabel(C,h,'Color','k','FontSize',25);
        else
            fprintf('%s: Rc=1 does not cross the interior of the grid.\n', ...
                names{panel});
        end

        xline(ax,annual_baselines(panel),':','Color',[1 1 1],'LineWidth',2);
        yline(ax,p_base,':','Color',[1 1 1],'LineWidth',2);
        set(ax,'FontSize',25,'TickLabelInterpreter','latex', ...
            'Box','off','TickDir','out','LineWidth',0.8);
        ax.Toolbar.Visible = 'off';
        xlabel(ax,xlabels{panel},'Interpreter','latex','FontSize',25);
        ylabel(ax,{'Annual maternal immunisation'; 'coverage, $p$'}, ...
            'Interpreter','latex','FontSize',25);
        ylim(ax,[0 1]);
        yticks(ax,0:0.2:1);
        yticklabels(ax,{'0\%','20\%','40\%','60\%','80\%','100\%'});
        if panel == 1
            xlim(ax,[0 1]);
            xticks(ax,0:0.2:1);
            xticklabels(ax,{'0\%','20\%','40\%','60\%','80\%','100\%'});
        else
            xlim(ax,[0 0.05]);
            xticks(ax,0:0.01:0.05);
            xticklabels(ax,{'0\%','1\%','2\%','3\%','4\%','5\%'});
        end
        % Force identical horizontal tick labels in both panels.
        % Explicit angles prevent MATLAB from rotating labels automatically.
        xtickangle(ax,0);
        ytickangle(ax,0);
        grid(ax,'on');
        hold(ax,'off');
        drawnow;

        % Export the whole fixed-size canvas, preserving matching margins.
        set(fig,'PaperUnits','inches','PaperPosition',[0 0 8.2 7.6], ...
            'PaperSize',[8.2 7.6],'PaperPositionMode','manual','InvertHardcopy','off');
        print(fig,fullfile(outdir_eps,[names{panel} '.eps']),'-depsc','-painters');
        print(fig,fullfile(outdir_png,[names{panel} '.png']),'-dpng','-r300');

        results.figures(panel) = fig;
        results.(names{panel}) = struct('psi_ann_grid',psi_ann_grid, ...
            'RCmat',RCmat,'baseline_error',baseline_error);
    end

    fprintf('EPS figures saved in: %s\n',outdir_eps);
    fprintf('PNG figures saved in: %s\n',outdir_png);
end
