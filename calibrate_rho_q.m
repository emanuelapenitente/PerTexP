function [rho_trans, q_trans, beta11, beta12, beta21, beta22, ...
          CI1, CI2, err1_rel, err2_rel, slack1, slack2, ...
          Jmin, R0_out, Rc_out, info] = ...
    calibrate_rho_q( ...
        par_base, contact_rates, x0, target1, target2, ...
        rho_bounds, q_trans_bounds, weights)
%CALIBRATE_RHO_Q Direct estimation of rho_trans and q_trans.
% C11,C12,C22 are WEEKLY Mistry rates. C21 is recomputed as
% (par_base.N1/par_base.N2)*C12, preserving contacts per infant.
% beta11=q*C11; beta12=q*C12; beta21=q*C21; beta22=rho*q*C22.
% rho_trans modifies transmission effectiveness only for non-infant contacts
% with other non-infants (beta22). It accounts for contact heterogeneity,
% not a biological reduction in susceptibility or infectivity.
% It is distinct from the reporting fractions rho1 and rho2.
% Defaults: rho_bounds=[0.01 0.999], q_trans_bounds=[1e-8 1], weights=[1 1].
% Requires MATLAB Optimization Toolbox (fmincon).
% Stage 1: multistart fmincon (SQP) with explicit parameter bounds,
% matching both Year-1 targets within relative tolerance 1e-6.
% Stage 2: minimize weighted squared relative errors subject to CI>=target,
% using a feasible grid followed by fmincon (SQP).
% Inputs and outputs are unchanged; x = [rho_trans; q_trans].
% An unsuccessful grid search does not prove mathematical infeasibility.

    if nargin < 6 || isempty(rho_bounds)
        rho_bounds = [0.01, 0.999];
    end
    if nargin < 7 || isempty(q_trans_bounds)
        q_trans_bounds = [1e-8, 1];
    end
    if nargin < 8 || isempty(weights)
        weights = [1, 1];
    end

    if exist('fmincon','file') ~= 2 || ~license('test','Optimization_Toolbox')
        error('calibrate_rho_q:MissingToolbox', ...
            'Calibration requires MATLAB Optimization Toolbox (fmincon).');
    end

    validate_bounds(q_trans_bounds, 'q_trans_bounds');
    if q_trans_bounds(2) > 1
        error('q_trans bounds must lie in (0,1].');
    end
    validateattributes(par_base.N1, {'numeric'}, {'scalar','finite','positive'});
    validateattributes(par_base.N2, {'numeric'}, {'scalar','finite','positive'});

    required_contact_fields = {'C11_week','C12_week','C22_week'};
    for kfield = 1:numel(required_contact_fields)
        if ~isfield(contact_rates, required_contact_fields{kfield})
            error('contact_rates is missing field %s.', required_contact_fields{kfield});
        end
    end

    C11 = contact_rates.C11_week;
    C12 = contact_rates.C12_week;
    C21 = (par_base.N1/par_base.N2)*C12;
    C22 = contact_rates.C22_week;

    validate_bounds(rho_bounds, 'rho_bounds');
    validate_bounds(q_trans_bounds, 'q_trans_bounds');

    if rho_bounds(2) > 1
        error('rho_trans is defined as a reduction factor and must satisfy rho <= 1.');
    end
    if any(~isfinite([C11,C12,C21,C22])) || any([C11,C12,C21,C22] <= 0)
        error('C11, C12, C21 and C22 must all be positive weekly contact rates.');
    end
    if numel(weights) ~= 2 || any(weights < 0) || all(weights == 0)
        error('weights must contain two nonnegative values, not both zero.');
    end

    % ===============================================================
    % STAGE 1: try to match both targets exactly.
    % ===============================================================
    rho_starts  = unique([0.05 0.10 0.20 0.35 0.50 0.65 0.80 0.95]);
    q_starts = unique([0.0001 0.001 0.005 0.01 0.03 0.1 0.3 0.8]);

    rho_starts = rho_starts(rho_starts > rho_bounds(1) & ...
                            rho_starts < rho_bounds(2));
    q_starts = q_starts(q_starts > q_trans_bounds(1) & ...
                              q_starts < q_trans_bounds(2));

    if isempty(rho_starts)
        rho_starts = mean(rho_bounds);
    end
    if isempty(q_starts)
        q_starts = mean(q_trans_bounds);
    end

    lb = [rho_bounds(1); q_trans_bounds(1)];
    ub = [rho_bounds(2); q_trans_bounds(2)];
    objfun = @(x) objective_physical( ...
        x, par_base, C11, C12, C21, C22, ...
        x0, target1, target2, weights);

    opts_exact = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'Display','off', ...
        'FiniteDifferenceType','central', ...
        'ConstraintTolerance',1e-10, ...
        'OptimalityTolerance',1e-12, ...
        'StepTolerance',1e-12, ...
        'MaxIterations',1000, ...
        'MaxFunctionEvaluations',4000);

    best_exact.J = Inf;
    best_exact.exitflag = NaN;
    best_exact.rho = NaN;
    best_exact.q_trans = NaN;
    best_exact.CI1 = NaN;
    best_exact.CI2 = NaN;

    nstarts = 0;

    for ir = 1:numel(rho_starts)
        for ib = 1:numel(q_starts)
            nstarts = nstarts + 1;

            xinit = [rho_starts(ir); q_starts(ib)];
            [xopt, ~, exitflag] = fmincon( ...
                objfun, xinit, [], [], [], [], lb, ub, [], opts_exact);
            rho_try = xopt(1);
            q_trans_try = xopt(2);

            [CI1_try, CI2_try, J_try] = metrics( ...
                rho_try, q_trans_try, par_base, ...
                C11, C12, C21, C22, x0, target1, target2, weights);

            if isfinite(J_try) && all(isfinite([CI1_try CI2_try])) && ...
               all(xopt >= lb) && all(xopt <= ub) && J_try < best_exact.J
                best_exact.J = J_try;
                best_exact.exitflag = exitflag;
                best_exact.rho = rho_try;
                best_exact.q_trans = q_trans_try;
                best_exact.CI1 = CI1_try;
                best_exact.CI2 = CI2_try;
            end
        end
    end

    exact_e1 = (best_exact.CI1-target1)/target1;
    exact_e2 = (best_exact.CI2-target2)/target2;

    exact_tol = 1e-6;  % max relative error tolerated as numerical exact match
    exact_found = all(isfinite([exact_e1 exact_e2])) && ...
        max(abs([exact_e1 exact_e2])) <= exact_tol;
    stage2_exitflag = [];
    stage2_refinement_accepted = false;

    if exact_found
        rho_trans = best_exact.rho;
        q_trans = best_exact.q_trans;
        CI1 = best_exact.CI1;
        CI2 = best_exact.CI2;
        Jmin = best_exact.J;

        calibration_mode = 'exact two-target match';
        fallback_method = 'not needed';

    else
        % ===========================================================
        % STAGE 2: constrained closest fit, always from above.
        % ===========================================================
        [rho0, q0, J0, CI10, CI20, grid_info] = ...
            best_feasible_grid_point( ...
                par_base, C11, C12, C21, C22, x0, target1, target2, ...
                rho_bounds, q_trans_bounds, weights);

        if isnan(rho0)
            error(['No feasible (rho_trans,q_trans) pair was found on the grid. ', ...
                   'Check bounds, targets, or grid resolution.']);
        end

        rho_trans = rho0;
        q_trans = q0;
        CI1 = CI10;
        CI2 = CI20;
        Jmin = J0;

        xinit = [rho0; q0];
        nonlcon = @(x) lower_bound_constraints( ...
            x, par_base, C11, C12, C21, C22, ...
            x0, target1, target2);

        opts_con = optimoptions('fmincon', ...
            'Algorithm','sqp', ...
            'Display','off', ...
            'FiniteDifferenceType','central', ...
            'ConstraintTolerance',1e-10, ...
            'OptimalityTolerance',1e-10, ...
            'StepTolerance',1e-12, ...
            'MaxIterations',500, ...
            'MaxFunctionEvaluations',5000);

        [xopt, ~, stage2_exitflag] = fmincon( ...
            objfun, xinit, [], [], [], [], lb, ub, nonlcon, opts_con);

        [CI1_try, CI2_try, J_try] = metrics( ...
            xopt(1), xopt(2), par_base, ...
            C11, C12, C21, C22, x0, target1, target2, weights);

        % Keep the feasible grid point unless refinement is feasible
        % and at least as good. Acceptance is checked independently
        % of the solver exit flag, which is reported in info.
        if all(isfinite([CI1_try CI2_try J_try])) && ...
           all(xopt >= lb) && all(xopt <= ub) && ...
           CI1_try >= target1 && CI2_try >= target2 && J_try <= Jmin
            rho_trans = xopt(1);
            q_trans = xopt(2);
            CI1 = CI1_try;
            CI2 = CI2_try;
            Jmin = J_try;
            stage2_refinement_accepted = true;
        end

        fallback_method = 'fmincon constrained optimization';

        calibration_mode = 'constrained closest fit from above';
    end

    % ===============================================================
    % Final reconstruction and diagnostics.
    % ===============================================================
    [beta11, beta12, beta21, beta22] = reconstruct_betas( ...
        rho_trans, q_trans, par_base, C11, C12, C21, C22);

    err1_rel = (CI1-target1)/target1;
    err2_rel = (CI2-target2)/target2;
    slack1 = CI1-target1;
    slack2 = CI2-target2;

    % If exact matching failed, the fallback solution MUST be from above.
    if ~exact_found
        tol1 = 1e-8*max(1,target1);
        tol2 = 1e-8*max(1,target2);

        if slack1 < -tol1 || slack2 < -tol2
            error(['Internal error: constrained fallback returned a ', ...
                   'solution below a required incidence lower bound.']);
        end
    end

    par = par_base;
    par.beta11 = beta11;
    par.beta12 = beta12;
    par.beta21 = beta21;
    par.beta22 = beta22;

    [~,~,~,~,R0_out,Rc_out] = RC_calc(par);

    info = struct();
    info.mode = calibration_mode;
    info.solver = 'fmincon';
    info.algorithm = 'sqp';
    info.best_exact_candidate_exitflag = best_exact.exitflag;
    info.stage2_exitflag = stage2_exitflag;
    info.stage2_refinement_accepted = stage2_refinement_accepted;
    info.fallback_method = fallback_method;
    info.exact_match_found = exact_found;
    info.exact_match_relative_tolerance = exact_tol;
    info.n_exact_multistarts = nstarts;
    info.best_exact_candidate_J = best_exact.J;
    info.best_exact_candidate_rho = best_exact.rho;
    info.best_exact_candidate_q_trans = best_exact.q_trans;
    info.rho_bounds = rho_bounds;
    info.q_trans_bounds = q_trans_bounds;
    info.weights = weights;
    info.free_parameters = {'rho_trans','q_trans'};
    info.C21_week_used = C21;
    info.contact_reciprocity_residual = par_base.N1*C12-par_base.N2*C21;

    if exist('grid_info','var')
        info.grid = grid_info;
    end
end


% =====================================================================
% Local helpers
% =====================================================================

function validate_bounds(b, name)
    if numel(b) ~= 2 || ~all(isfinite(b)) || b(1) <= 0 || b(2) <= b(1)
        error('%s must be [lower upper] with 0 < lower < upper.', name);
    end
end


function [CI1, CI2, J] = metrics( ...
    rho_trans, q_trans, par_base, ...
    C11, C12, C21, C22, x0, target1, target2, weights)

    [beta11, beta12, beta21, beta22] = reconstruct_betas( ...
        rho_trans, q_trans, par_base, C11, C12, C21, C22);

    par = par_base;
    par.beta11 = beta11;
    par.beta12 = beta12;
    par.beta21 = beta21;
    par.beta22 = beta22;

    [~,~,~,inc1,inc2,~,~,~,~,~] = discrete_system(0,52,x0,par);

    CI1 = sum(inc1);
    CI2 = sum(inc2);

    e1 = (CI1-target1)/target1;
    e2 = (CI2-target2)/target2;

    J = weights(1)*e1^2 + weights(2)*e2^2;
end


function J = objective_physical( ...
    x, par_base, C11, C12, C21, C22, ...
    x0, target1, target2, weights)

    [~,~,J] = metrics( ...
        x(1), x(2), par_base, C11, C12, C21, C22, ...
        x0, target1, target2, weights);
end


function [c,ceq] = lower_bound_constraints( ...
    x, par_base, C11, C12, C21, C22, x0, target1, target2)

    [CI1,CI2,~] = metrics( ...
        x(1), x(2), par_base, C11, C12, C21, C22, ...
        x0, target1, target2, [1 1]);

    % fmincon uses c(x) <= 0.
    c = [target1-CI1; target2-CI2];
    ceq = [];
end


function [rho_best, q_best, J_best, CI1_best, CI2_best, info] = ...
    best_feasible_grid_point( ...
        par_base, C11, C12, C21, C22, x0, target1, target2, ...
        rho_bounds, q_trans_bounds, weights)

    % Moderate grid: used only if exact matching was not found.
    nrho = 90;
    nq = 180;

    rho_grid = linspace(rho_bounds(1), rho_bounds(2), nrho);
    q_grid = logspace(log10(q_trans_bounds(1)), log10(q_trans_bounds(2)), nq);

    rho_best = NaN;
    q_best = NaN;
    J_best = Inf;
    CI1_best = NaN;
    CI2_best = NaN;
    nfeasible = 0;

    for ir = 1:nrho
        rho_try = rho_grid(ir);

        for ib = 1:nq
            q_try = q_grid(ib);

            [CI1_try, CI2_try, J_try] = metrics( ...
                rho_try, q_try, par_base, ...
                C11, C12, C21, C22, x0, target1, target2, weights);

            if CI1_try >= target1 && CI2_try >= target2
                nfeasible = nfeasible + 1;

                if J_try < J_best
                    rho_best = rho_try;
                    q_best = q_try;
                    J_best = J_try;
                    CI1_best = CI1_try;
                    CI2_best = CI2_try;
                end
            end
        end
    end

    info = struct();
    info.nrho = nrho;
    info.nq = nq;
    info.nfeasible = nfeasible;
end


function [beta11, beta12, beta21, beta22] = reconstruct_betas( ...
    rho_trans, q_trans, par_base, C11, C12, C21, C22)
    beta11 = q_trans*C11;
    beta12 = q_trans*C12;
    beta21 = q_trans*C21;
    beta22 = rho_trans*q_trans*C22;
end
